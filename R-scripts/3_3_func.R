library(caret)
library(dplyr)
library(terra)

# u_graph <- function(x, y, type) {
# 
#   plot(x, y, type = "l", ylim = c(0, 1), col = "blue", ylab = "Condition")
#   
#   
# }
### pkg/R/goofcat.R
# Purpose        : Goodness of fit statistics for categorical varibles
# Maintainer     : Brendan Malone (brendan.malone@sydney.edu.au); 
# Contributions  : 
# Status         : working
# Note           : 

clean_temp <- function(){file.remove(list.files(tempdir(), full.names = T, pattern = ".tif"))}

pre.ab <- function(x) {
  x <- ifelse(x > 0.5, 1, 0)
}

inAustralia <- function(soilDF){
  if(!is.data.frame(soilDF)){return(NULL)}
  bboxExt <- extent(110,155,-45,-9)
  idxs <- which(soilDF$Longitude >= bboxExt@xmin & soilDF$Longitude <= bboxExt@xmax & soilDF$Latitude >= bboxExt@ymin & soilDF$Latitude <= bboxExt@ymax)
  outdf <- soilDF[idxs, ]
}

removeDuplicateRecords <- function(soilDF){
  noDups <- soilDF[!duplicated(soilDF[,1:15]),] ## I don't include the extraction time and last variables
  return(noDups)
}

func.trans.grav.vol <- function(x,bd){
  CFvol <- ifelse(test= (is.na(bd)|is.na(x)), yes=NA, no= x*(bd/2.65))
  return(CFvol)
}

goofcat<- function(observed = NULL, predicted = NULL, conf.mat, imp=FALSE){
  
  if (imp==TRUE){
    if(class(conf.mat)!="matrix"){
      stop("Entered data is NOT a matrix")}
    if(nrow(conf.mat)!= ncol(conf.mat)) {
      stop("Entered data is NOT a confusion matrix")}
    
    else {OA<- ceiling(sum(diag(conf.mat))/sum(colSums(conf.mat)) * 100)
    PA<- ceiling(diag(conf.mat)/colSums(conf.mat) * 100)
    UA<- ceiling(diag(conf.mat)/rowSums(conf.mat) * 100)
    
    PE_mat <- matrix(NA, ncol = 1, nrow = length(rowSums(conf.mat)))
    for (i in 1:length(rowSums(conf.mat))) {
      PE_mat[i, 1] <- (rowSums(conf.mat)[i]/sum(colSums(conf.mat))) * (colSums(conf.mat)[i]/sum(colSums(conf.mat)))}
    KS <- (sum(diag(conf.mat))/sum(colSums(conf.mat)) - sum(PE_mat))/(1 - sum(PE_mat))}}
  
  if (imp==FALSE) {obsMat<- table(observed,observed)
  df<- data.frame(observed, predicted)
  names(df)<- c("observed", "predicted")
  #make a confusion matrix
  cfuM<- function(df,obsMat){
    c.Mat<- as.matrix(obsMat)
    snames1<- c(colnames(c.Mat))
    for (i in 1:nrow(c.Mat)){
      for (j in 1:nrow(c.Mat)){
        c.Mat[j,i]<- nrow(subset(df, df$observed ==snames1[i]  & df$predicted ==snames1[j]))}}
    fmat<- matrix(NA, nrow=nrow(c.Mat), ncol=ncol(c.Mat))
    rownames(fmat)<- rownames(c.Mat)
    colnames(fmat)<- colnames(c.Mat)
    for (i in 1:nrow(c.Mat)){
      fmat[i,]<- c(c.Mat[i,])}
    return(fmat)}
  conf.mat<- cfuM(df, obsMat)
  
  OA<- ceiling(sum(diag(conf.mat))/sum(colSums(conf.mat)) * 100)
  PA<- ceiling(diag(conf.mat)/colSums(conf.mat) * 100)
  UA<- ceiling(diag(conf.mat)/rowSums(conf.mat) * 100)
  
  PE_mat <- matrix(NA, ncol = 1, nrow = length(rowSums(conf.mat)))
  for (i in 1:length(rowSums(conf.mat))) {
    PE_mat[i, 1] <- (rowSums(conf.mat)[i]/sum(colSums(conf.mat))) * (colSums(conf.mat)[i]/sum(colSums(conf.mat)))}
  KS <- (sum(diag(conf.mat))/sum(colSums(conf.mat)) - sum(PE_mat))/(1 - sum(PE_mat))}
  retval<- list(conf.mat,OA, PA, UA, KS)
  names(retval)<- c("confusion_matrix", "overall_accuracy", "producers_accuracy", "users_accuracy", "kappa")
  return(retval)}


lm.utg <- function(r, x, y,mode) {
  if (length(x)>3){
    x = x[2:3]
    y = y[2:3]
  } else if(length(x)==3){
    if (mode == "max"){
      x = x[-length(x)]
      y = y[-length(y)]
    } else if (mode=='min'){
      x = x[-1]
      y = y[-1]
    }
  }
  
  # create the linear regression model to make the values
  dat <- data.frame(cbind(x, y))
  model <- lm(y ~ ., data = dat)
  
  r.utg <- app(r, fun = function(cc) {
    cc <- model$coefficients[2] * cc + model$coefficients[1]
    return(cc)
  })
  
  if(mode == "max"){
    r.utg[r > max(x)] <- 1
    r.utg[r < min(x)] <- 0
  } else if(mode=='min') {
    r.utg[r < min(x)] <- 1
    r.utg[r > max(x)] <- 0
  }
  
  gc()
  return(r.utg)
}

mod.utg <- function(r, x, y, max = T) {
  # create the linear regression model to make the values
  dat <- data.frame(cbind(x, y))
  model <- lm(y ~ ., data = dat)
  
  # get the raster values, and predicted values from linear model
  df <- data.frame(values(r))
  df$val <- predict(model, newdata = data.frame(x = unlist(df)))
  
  # create utility graph to predict on the raster
  model <- lm(val ~ ., data = df)
  r.utg <- predict(r, model)
  if (max == T) {
    r.utg[r > max(x)] <- 1
  } else {
    r.utg[r < min(x)] <- 1
  }
  # plot(r.utg)
  gc()
  return(r.utg)
}

# Get the existing method for the ranger() function
ranger_type = getModelInfo(model = "ranger")
ranger_type

# Change the parameters that may be tuned to include num.trees
ranger_type$ranger$parameters = data.frame("parameter" = c("mtry", "num.trees","splitrule","min.node.size"),
                                           "class" = c("numeric", "numeric"),
                                           "label" = c("Selected Predictors", 
                                                       "Number of Trees"))

# Edit the model fit function to include a num.trees argument in the call to the ranger function()
ranger_type$ranger$fit = function (x, y, wts, param, lev, last, classProbs, 
                                   ...) {
  if ((!is.data.frame(x)) || dplyr::is.tbl(x)) 
    x <- as.data.frame(x)
  x$.outcome <- y
  if (!is.null(wts)) {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = param$mtry, num.trees = param$num.trees,
                          probability = classProbs, case.weights = wts, ...)
  }
  else {
    out <- ranger::ranger(dependent.variable.name = ".outcome", 
                          data = x, mtry = param$mtry, num.trees = param$num.trees, 
                          write.forest = TRUE, probability = classProbs, ...)
  }
  if (!last) 
    out$y <- y
  out
}

goof<-function (observed, predicted, coefficient = c("R2", "concordance", 
                                                     "MSE", "RMSE", "bias", "RPIQ"), plot = TRUE, ...) 
{
  if (any(!coefficient %in% c("R2", "concordance", 
                              "MSE", "RMSE", "bias", "MSEc", 
                              "RMSEc", "RPD", "RPIQ"))) 
    stop("Please choose a valid coefficient")
  rLM <- lm(predicted ~ observed)
  R2 <- as.matrix(summary(rLM)$adj.r.squared)
  SEP2 <- mean((observed - predicted)^2)
  SEP <- sqrt(SEP2)
  bias <- mean(predicted) - mean(observed)
  SEP2c <- sum(((predicted - bias - observed)^2)/length(observed))
  SEPc <- sqrt(SEP2c)
  RPD <- sd(observed)/SEP
  IQ <- c(quantile(observed))[4] - c(quantile(observed))[2]
  RPIQ <- IQ/SEP
  mx <- mean(observed)
  my <- mean(predicted)
  s2x <- var(observed)
  s2y <- var(predicted)
  sxy <- mean((observed - mx) * (predicted - my))
  ccc <- 2 * sxy/(s2x + s2y + (mx - my)^2)
  if (plot) {
    plot(observed, predicted, ylim = c(min(c(observed, predicted)), 
                                       max(c(observed, predicted))), xlim = c(min(c(observed, 
                                                                                    predicted)), max(c(observed, predicted))), asp = 1, 
         ...)
    abline(a = 0, b = 1, col = "brown4")
  }
  coefs_tmp <- data.frame(R2 = R2, concordance = ccc, MSE = SEP2, 
                          RMSE = SEP, bias = bias, RPIQ = RPIQ, row.names = NULL)
  gf <- data.frame(coefs_tmp[, coefficient])
  
  
  gf
}


# Remove observations less than cutoffs value -----------------------------
na_fill <- function(df, cutoffs, col_names){
  for (cc in col_names){
    
    # Subset the original data frame with the targeted columns
    df1 = df[,cc]
    
    #Apply the replacement with NA in all columns in the 2nd data frame.
    if (length(cutoffs)>1){
    df1[df1 <= cutoffs[1]] = NA
    df1[df1 >= cutoffs[2]] = NA
    } else{
    df1[df1 >= cutoffs] = NA}
    # Assign the resulted column in the 2nd data frame to the original data frame
    df[,cc] = df1
  }
  return(df)
}

ea_spline <- function (obj, var.name, lam = 0.1, d = c(0, 5, 15, 30, 60, 100, 
                                                       200), vlow = 0, vhigh = 1000, show.progress = TRUE) 
{
  if (is(obj, "SoilProfileCollection") == TRUE) {
    depthcols = obj@depthcols
    idcol = obj@idcol
    obj@horizons = obj@horizons[, c(idcol, depthcols, var.name)]
    d.dat <- as.data.frame(obj@horizons)
  }
  if (is(obj, "data.frame") == TRUE) {
    d.dat <- as.data.frame(obj[, c(1:3, which(colnames(obj) == 
                                                var.name))])
  }
  if (is(obj, "data.frame") == FALSE & is(obj, "SoilProfileCollection") == 
      FALSE) {
    stop("ERROR:: Data must be either class data.frame or SoilProfileCollection")
  }
  mxd <- max(d)
  sp_dat <- split(d.dat, d.dat[, 1])
  m_fyfit <- matrix(NA, ncol = length(c(1:mxd)), nrow = length(sp_dat))
  yave <- matrix(NA, ncol = length(d), nrow = length(sp_dat))
  sse <- matrix(NA, ncol = length(lam), nrow = 1)
  sset <- matrix(NA, ncol = 2, nrow = length(sp_dat))
  mat_id <- d.dat[0, ]
  dave <- d.dat[1, ]
  dave$predicted <- 0
  dave$FID <- 0
  dave <- dave[0, ]
  if (show.progress) 
    pb <- txtProgressBar(min = 0, max = length(sp_dat), style = 3)
  cnt <- 1
  for (st in 1:length(sp_dat)) {
    subs <- sp_dat[[st]]
    subs <- as.matrix(subs)
    mat_id[st, 1] <- subs[1, 1]
    ir <- c(1:nrow(subs))
    ir <- as.matrix(t(ir))
    u <- as.numeric(subs[ir, 2])
    u <- as.matrix(t(u))
    v <- as.numeric(subs[ir, 3])
    v <- as.matrix(t(v))
    y <- as.numeric(subs[ir, 4])
    y <- as.matrix(t(y))
    n <- length(y)
    d <- t(d)
    if (n == 1) {
      dave[cnt:(cnt - 1 + nrow(subs)), 1:4] <- subs
      dave[cnt:(cnt - 1 + nrow(subs)), 5] <- y
      dave[cnt:(cnt - 1 + nrow(subs)), 6] <- st
      xfit <- as.matrix(t(c(1:mxd)))
      nj <- max(v)
      if (nj > mxd) {
        nj <- mxd
      }
      yfit <- xfit
      yfit[, 1:nj] <- y
      if (nj < mxd) {
        yfit[, (nj + 1):mxd] = -9999
      }
      m_fyfit[st, ] <- yfit
      nd <- length(d) - 1
      dl <- d + 1
      for (cj in 1:nd) {
        xd1 <- dl[cj]
        xd2 <- dl[cj + 1] - 1
        if (nj >= xd1 & nj <= xd2) {
          xd2 <- nj - 1
          yave[st, cj] <- mean(yfit[, xd1:xd2])
        }
        else {
          yave[st, cj] <- mean(yfit[, xd1:xd2])
        }
        yave[st, cj + 1] <- max(v)
      }
      cnt <- cnt + nrow(subs)
      sset[st, 1:2] <- 0
    }
    else {
      dave[cnt:(cnt - 1 + nrow(subs)), 1:4] <- subs
      dave[cnt:(cnt - 1 + nrow(subs)), 6] <- st
      np1 <- n + 1
      nm1 <- n - 1
      delta <- v - u
      del <- c(u[2:n], u[n]) - v
      r <- matrix(0, ncol = nm1, nrow = nm1)
      for (dig in 1:nm1) {
        r[dig, dig] <- 1
      }
      for (udig in 1:nm1 - 1) {
        r[udig, udig + 1] <- 1
      }
      d2 <- matrix(0, ncol = nm1, nrow = nm1)
      diag(d2) <- delta[2:n]
      r <- d2 %*% r
      r <- r + t(r)
      d1 <- matrix(0, ncol = nm1, nrow = nm1)
      diag(d1) <- delta[1:nm1]
      d3 <- matrix(0, ncol = nm1, nrow = nm1)
      diag(d3) <- del[1:nm1]
      r <- r + 2 * d1 + 6 * d3
      q <- matrix(0, ncol = n, nrow = n)
      for (dig in 1:n) {
        q[dig, dig] <- -1
      }
      for (udig in 1:n - 1) {
        q[udig, udig + 1] <- 1
      }
      q <- q[1:nm1, 1:n]
      dim.mat <- matrix(q[], ncol = length(1:n), nrow = length(1:nm1))
      rinv <- try(solve(r), TRUE)
      if (is.matrix(rinv)) {
        ind <- diag(n)
        pr.mat <- matrix(0, ncol = length(1:nm1), nrow = length(1:n))
        pr.mat[] <- 6 * n * lam
        fdub <- pr.mat * t(dim.mat) %*% rinv
        z <- fdub %*% dim.mat + ind
        sbar <- solve(z, t(y))
        dave[cnt:(cnt - 1 + nrow(subs)), 5] <- sbar
        cnt <- cnt + nrow(subs)
        b <- 6 * rinv %*% dim.mat %*% sbar
        b0 <- rbind(0, b)
        b1 <- rbind(b, 0)
        gamma <- (b1 - b0)/t(2 * delta)
        alfa <- sbar - b0 * t(delta)/2 - gamma * t(delta)^2/3
        xfit <- as.matrix(t(c(1:mxd)))
        nj <- max(v)
        if (nj > mxd) {
          nj <- mxd
        }
        yfit <- xfit
        for (k in 1:nj) {
          xd <- xfit[k]
          if (xd < u[1]) {
            p <- alfa[1]
          }
          else {
            for (its in 1:n) {
              if (its < n) {
                tf2 = as.numeric(xd > v[its] & xd < u[its + 
                                                        1])
              }
              else {
                tf2 <- 0
              }
              if (xd >= u[its] & xd <= v[its]) {
                p = alfa[its] + b0[its] * (xd - u[its]) + 
                  gamma[its] * (xd - u[its])^2
              }
              else if (tf2) {
                phi = alfa[its + 1] - b1[its] * (u[its + 
                                                     1] - v[its])
                p = phi + b1[its] * (xd - v[its])
              }
            }
          }
          yfit[k] = p
        }
        if (nj < mxd) {
          yfit[, (nj + 1):mxd] = NA
        }
        yfit[which(yfit > vhigh)] <- vhigh
        yfit[which(yfit < vlow)] <- vlow
        m_fyfit[st, ] <- yfit
        nd <- length(d) - 1
        dl <- d + 1
        for (cj in 1:nd) {
          xd1 <- dl[cj]
          xd2 <- dl[cj + 1] - 1
          if (nj >= xd1 & nj <= xd2) {
            xd2 <- nj - 1
            yave[st, cj] <- mean(yfit[, xd1:xd2])
          }
          else {
            yave[st, cj] <- mean(yfit[, xd1:xd2])
          }
          yave[st, cj + 1] <- max(v)
        }
        rmse <- sqrt(sum((t(y) - sbar)^2)/n)
        rmseiqr <- rmse/IQR(y)
        sset[st, 1] <- rmse
        sset[st, 2] <- rmseiqr
      }
    }
    if (show.progress) {
      setTxtProgressBar(pb, st)
    }
  }
  if (show.progress) {
    close(pb)
  }
  yave <- as.data.frame(yave)
  jmat <- matrix(NA, ncol = 1, nrow = length(d))
  for (i in 1:length(d) - 1) {
    a1 <- paste(d[i], d[i + 1], sep = "-")
    a1 <- paste(a1, "cm", sep = " ")
    jmat[i] <- a1
  }
  jmat[length(d)] <- "soil depth"
  for (jj in 1:length(jmat)) {
    names(yave)[jj] <- jmat[jj]
  }
  yave <- cbind(mat_id[, 1], yave)
  names(yave)[1] <- "id"
  sset <- as.data.frame(sset)
  names(sset) <- c("rmse", "rmseiqr")
  retval <- list(harmonised = yave, obs.preds = dave, splineFitError = sset, 
                 var.1cm = t(m_fyfit))
  return(retval)
}

get_sd <-function(ras,ras95){
  ras.m = global(ras, "mean", na.rm = T)$mean
  temp = (ras95-ras.m)/(1.645*ras)
  temp[temp<0] =0
  return(temp)
}

clean_up<-function(){
  .rs.restartR()
  folders <- dir(Sys.getenv("TEMP"), pattern = "Rtmp", full.names = TRUE)
  unlink(folders, recursive = TRUE, force = TRUE, expand = TRUE)}
