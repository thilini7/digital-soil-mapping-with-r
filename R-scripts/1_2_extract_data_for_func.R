#' ea_spline: Fits a mass-preserving spline to soil profile data
#'
#' Fits equal-area smoothing splines to numeric values across soil depth intervals (horizons)
#' for each profile, preserving the original mass (integral) of observations.
#' Accepts either a data.frame or SoilProfileCollection.
#'
#' @param obj A `data.frame` or `SoilProfileCollection`. Must include profile ID, upper and lower depths, and the target variable.
#' @param var.name Character. The name of the target numeric variable to be spline-smoothed.
#' @param lam Numeric. Spline smoothing parameter (default = 0.1).
#' @param d Numeric vector. Target standard depth intervals (default = c(0,5,15,30,60,100,200)).
#' @param vlow Numeric. Lower bound for truncating predictions (default = 0).
#' @param vhigh Numeric. Upper bound for truncating predictions (default = 1000).
#' @param show.progress Logical. Whether to show a progress bar (default = TRUE).
#'
#' @return A list with:
#' \item{harmonised}{Data frame of spline-estimated values at standard depths.}
#' \item{obs.preds}{Data frame of original and predicted values.}
#' \item{splineFitError}{RMSE and RMSE/IQR for each profile.}
#' \item{var.1cm}{Matrix of 1cm spline predictions for each profile.}
#'
ea_spline <- function(obj, var.name, lam = 0.1,
                      d = c(0, 5, 15, 30, 60, 100, 200),
                      vlow = 0, vhigh = 1000,
                      show.progress = TRUE) {
  # --- Validate and prepare input ---
  if (inherits(obj, "SoilProfileCollection")) {
    depthcols <- obj@depthcols
    idcol <- obj@idcol
    d.dat <- as.data.frame(obj@horizons[, c(idcol, depthcols, var.name)])
  } else if (is.data.frame(obj)) {
    d.dat <- as.data.frame(obj[, c(1:3, which(colnames(obj) == var.name))])
  } else {
    stop("Input must be a data.frame or SoilProfileCollection.")
  }
  
  max.depth <- max(d)
  profile.list <- split(d.dat, d.dat[, 1])
  n.profiles <- length(profile.list)
  
  spline.1cm <- matrix(NA, nrow = n.profiles, ncol = max.depth)
  spline.avg <- matrix(NA, nrow = n.profiles, ncol = length(d))
  fit.errors <- matrix(NA, nrow = n.profiles, ncol = 2)
  
  mat_id <- d.dat[0, ]
  obs.pred <- d.dat[1, ]
  obs.pred$predicted <- 0
  obs.pred$FID <- 0
  obs.pred <- obs.pred[0, ]
  
  if (show.progress) pb <- txtProgressBar(min = 0, max = n.profiles, style = 3)
  
  for (i in seq_len(n.profiles)) {
    profile <- as.matrix(profile.list[[i]])
    n.horizons <- nrow(profile)
    mat_id[i, 1] <- profile[1, 1]
    
    u <- as.numeric(profile[, 2])
    v <- as.numeric(profile[, 3])
    y <- as.numeric(profile[, 4])
    
    # --- Handle single-horizon profile ---
    if (n.horizons == 1) {
      pred <- matrix(rep(y, max.depth), nrow = 1)
      # Safely fill NA beyond LowerDepth
      if (!is.na(v[1]) && v[1] < max.depth) {
        idx <- seq(v[1] + 1, max.depth)
        pred[1, idx] <- NA
      }
      spline.1cm[i, ] <- pred
      
      # Average over target depths
      for (j in seq_len(length(d) - 1)) {
        d.start <- d[j] + 1
        d.end <- d[j + 1]
        spline.avg[i, j] <- mean(pred[, d.start:d.end], na.rm = TRUE)
      }
      spline.avg[i, length(d)] <- max(v)
      fit.errors[i, ] <- 0
      
      obs.pred <- rbind(obs.pred, cbind(profile, predicted = y, FID = i))
      
    } else {
      # --- Fit spline ---
      delta <- v - u
      del <- c(u[-1], u[n.horizons]) - v
      
      nm1 <- n.horizons - 1
      
      # r matrix
      r <- diag(1, nm1)
      r[row(r) == col(r) - 1] <- 1
      # d2 matrix
      if (nm1 > 1) {
        d2 <- diag(delta[-1])
      } else {
        d2 <- matrix(delta[2], nrow = 1, ncol = 1)
      }
      
      # enforce matrix diag
      Ddelta <- diag(delta[1:nm1], nrow = nm1, ncol = nm1)
      Ddel   <- diag(del[1:nm1],   nrow = nm1, ncol = nm1)
      
      # update r
      r <- d2 %*% r + t(d2 %*% r) + 2 * Ddelta + 6 * Ddel
      
      # q matrix
      q <- diag(-1, n.horizons)
      q[row(q) == col(q) - 1] <- 1
      q <- q[1:nm1, , drop = FALSE]  # keep matrix
      
      # safe inverse
      rinv <- tryCatch(solve(r), error = function(e) NULL)
      if (is.null(rinv)) {
        cat("⚠️ Singular r, skipping profile", i, "\n")
        next
      }
      rinv <- as.matrix(rinv) # enforce matrix
      
      # check multiplication
      tmp <- tryCatch(t(q) %*% rinv %*% q, error = function(e) {
        cat("⚠️ Error in t(q) %*% rinv %*% q:", e$message, "\n")
        NULL
      })
      if (is.null(tmp)) next
      
      # build z
      z <- 6 * n.horizons * lam * tmp + diag(n.horizons)
      
      # fit spline parameters
      sbar <- solve(z, y)
      b <- 6 * rinv %*% q %*% sbar
      b0 <- c(0, b)
      b1 <- c(b, 0)
      gamma <- (b1 - b0) / (2 * delta)
      alpha <- sbar - b0 * delta / 2 - gamma * delta^2 / 3
      
      # Predict at 1 cm increments
      xfit <- seq_len(max.depth)
      pred <- rep(NA, max.depth)
      for (k in xfit) {
        xd <- k
        if (xd < u[1]) {
          pred[k] <- alpha[1]
        } else {
          for (j in seq_len(n.horizons)) {
            if (xd >= u[j] && xd <= v[j]) {
              pred[k] <- alpha[j] + b0[j] * (xd - u[j]) + gamma[j] * (xd - u[j])^2
              break
            } else if (j < n.horizons && xd > v[j] && xd < u[j + 1]) {
              phi <- alpha[j + 1] - b1[j] * (u[j + 1] - v[j])
              pred[k] <- phi + b1[j] * (xd - v[j])
              break
            }
          }
        }
      }
      
      pred[pred > vhigh] <- vhigh
      pred[pred < vlow]  <- vlow
      spline.1cm[i, ] <- pred
      
      for (j in seq_len(length(d) - 1)) {
        d.start <- d[j] + 1
        d.end <- d[j + 1]
        spline.avg[i, j] <- mean(pred[d.start:d.end], na.rm = TRUE)
      }
      spline.avg[i, length(d)] <- max(v)
      
      rmse <- sqrt(mean((y - sbar)^2))
      rmse_iqr <- rmse / IQR(y)
      fit.errors[i, ] <- c(rmse, rmse_iqr)
      
      obs.pred <- rbind(obs.pred, cbind(profile, predicted = sbar, FID = i))
    }
    
    if (show.progress) setTxtProgressBar(pb, i)
  }
  
  if (show.progress) close(pb)
  
  # --- Final formatting ---
  spline.avg <- as.data.frame(spline.avg)
  colnames(spline.avg) <- c(
    paste0(d[-length(d)], "-", d[-1], "cm"), "soil depth"
  )
  spline.avg <- cbind(id = mat_id[, 1], spline.avg)
  
  fit.errors <- as.data.frame(fit.errors)
  colnames(fit.errors) <- c("rmse", "rmseiqr")
  
  return(list(
    harmonised = spline.avg,
    obs.preds = obs.pred,
    splineFitError = fit.errors,
    var.1cm = t(spline.1cm)
  ))
}
