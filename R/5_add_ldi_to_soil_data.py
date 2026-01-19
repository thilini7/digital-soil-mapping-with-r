#!/usr/bin/env python3
"""
Script to add LDI (Land Disturbance Index) values to soil point data
by performing spatial join with NSW Land Use shapefile.

LDI Classification (from Gray et al., 2015):
- LDI 1: No effective disturbance (e.g., national park, nature reserve)
- LDI 2: Limited disturbance, minor native vegetation clearing (e.g., selective logging, production forestry)
- LDI 3: Moderate disturbance, moderate native vegetation clearing, light grazing in woodland, hardwood plantation
- LDI 4: High disturbance, complete native vegetation clearing (e.g., native and improved pasture, softwood plantation)
- LDI 5: Very high disturbance (e.g., improved pasture with moderate cropping, orchards, viticulture)
- LDI 6: Extreme disturbance, predominant cropping (rain-fed or irrigated)
"""

import os
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point

# Define paths
BASE_DIR = "/Users/neo/Development/Thilini-git/digital-soil-mapping-with-r"
SHAPEFILE_PATH = os.path.join(BASE_DIR, "Data/data_in/shp_files/NSW_ACT_Landuse_merged.shp")
CSV_INPUT_PATH = os.path.join(BASE_DIR, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined.csv")
CSV_OUTPUT_PATH = os.path.join(BASE_DIR, "Data/data_out/ansis_lab_measurements/All Sites/ANSIS_combined_with_LDI.csv")

def main():
    print("=" * 60)
    print("Adding LDI (Land Disturbance Index) to soil point data")
    print("=" * 60)
    
    # Step 1: Load the soil point data
    print("\n1. Loading soil point data...")
    soil_df = pd.read_csv(CSV_INPUT_PATH)
    print(f"   Loaded {len(soil_df)} records")
    print(f"   Columns: {list(soil_df.columns)}")
    
    # Check for coordinate columns
    if 'Longitude' not in soil_df.columns or 'Latitude' not in soil_df.columns:
        raise ValueError("CSV must have 'Longitude' and 'Latitude' columns")
    
    # Get unique locations to reduce processing
    print("\n2. Extracting unique locations...")
    unique_locs = soil_df[['Longitude', 'Latitude']].drop_duplicates()
    print(f"   Found {len(unique_locs)} unique locations")
    
    # Step 2: Convert to GeoDataFrame
    print("\n3. Converting to spatial points...")
    geometry = [Point(xy) for xy in zip(unique_locs['Longitude'], unique_locs['Latitude'])]
    points_gdf = gpd.GeoDataFrame(unique_locs, geometry=geometry, crs="EPSG:4283")  # GDA94
    
    # Step 3: Load the land use shapefile
    print("\n4. Loading NSW Land Use shapefile (this may take a while)...")
    landuse_gdf = gpd.read_file(SHAPEFILE_PATH)
    print(f"   Loaded {len(landuse_gdf)} land use polygons")
    print(f"   Shapefile CRS: {landuse_gdf.crs}")
    print(f"   Available columns: {list(landuse_gdf.columns)}")
    
    # Check LDI values in shapefile
    print("\n5. LDI value distribution in shapefile:")
    print(landuse_gdf['LDI'].value_counts().sort_index())
    
    # Ensure CRS match
    if points_gdf.crs != landuse_gdf.crs:
        print(f"\n   Reprojecting points from {points_gdf.crs} to {landuse_gdf.crs}")
        points_gdf = points_gdf.to_crs(landuse_gdf.crs)
    
    # Step 4: Spatial join
    print("\n6. Performing spatial join...")
    joined_gdf = gpd.sjoin(points_gdf, landuse_gdf[['geometry', 'LDI', 'Secondary', 'Tertiary']], 
                           how='left', predicate='within')
    
    # Handle duplicates from sjoin (point in multiple polygons)
    joined_gdf = joined_gdf.drop_duplicates(subset=['Longitude', 'Latitude'])
    
    print(f"   Joined {len(joined_gdf)} unique locations")
    print(f"   Points with LDI: {joined_gdf['LDI'].notna().sum()}")
    print(f"   Points without LDI (outside polygons): {joined_gdf['LDI'].isna().sum()}")
    
    # Step 5: Merge back to original data
    print("\n7. Merging LDI values back to original data...")
    ldi_lookup = joined_gdf[['Longitude', 'Latitude', 'LDI', 'Secondary', 'Tertiary']].copy()
    ldi_lookup = ldi_lookup.rename(columns={'Secondary': 'LandUse_Secondary', 'Tertiary': 'LandUse_Tertiary'})
    
    # Merge with original data
    result_df = soil_df.merge(ldi_lookup, on=['Longitude', 'Latitude'], how='left')
    
    # Step 6: Save output
    print("\n8. Saving results...")
    result_df.to_csv(CSV_OUTPUT_PATH, index=False)
    print(f"   Saved to: {CSV_OUTPUT_PATH}")
    
    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Total records: {len(result_df)}")
    print(f"Records with LDI: {result_df['LDI'].notna().sum()}")
    print(f"Records without LDI: {result_df['LDI'].isna().sum()}")
    print("\nLDI distribution in output:")
    print(result_df['LDI'].value_counts().sort_index())

if __name__ == "__main__":
    main()
