# -*- coding: utf-8 -*-
"""
Script to extract montly weather data from the Info-climat gridded daily data.

http://www.environnement.gouv.qc.ca/climat/surveillance/produits.htm
"""

import os.path as osp
import numpy as np
import netCDF4
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point
import datetime as dt

infoclimat_dirname = "D:/Data/MeteoGrilleDaily"
project_dirname = osp.dirname(__file__)

# %% Load study area shapefile

print("Loading shapefile of study area... ", end='')
shpfilename = osp.join(
    project_dirname, "data_extraction_limit", "data_extraction_limit.shp")
zone_gdf = gpd.read_file(shpfilename)
print('done')

# %% Create attribute table

print("Creating an attribute table... ", end='')

# Get grid coordinates from one of the netCDF file.
netcdf_dset = netCDF4.Dataset(
    osp.join(infoclimat_dirname, "GCQ_v2_2016.nc"), 'r+')

latitudes = np.array(netcdf_dset['lat'])
longitudes = np.array(netcdf_dset['lon'])

# Generate an attribute table with the coordinates of the info-climat grid.
data = []
geometry = []
for i, lat in enumerate(latitudes):
    for j, lon in enumerate(longitudes):
        geometry.append(Point((lon, lat)))
        data.append([i, j, lat, lon])
data = pd.DataFrame(
    data, columns=['lat_idx', 'lon_idx', 'latitude', 'longitude'])
crs = "+proj=longlat +ellps=GRS80 +datum=NAD83 +towgs84=0,0,0,0,0,0,0 +no_defs"
sta_gdf = gpd.GeoDataFrame(data, crs=crs, geometry=geometry)

# Change the coordinate reference system of the info-climat grid to that
# of the study area.
sta_gdf = sta_gdf.to_crs(zone_gdf.crs)

# Define the points of the grid that are in the study area.
sta_gdf['InZone'] = sta_gdf['geometry'].within(
    zone_gdf.loc[0]['geometry']).astype(int)

print('done')

# %% Extract daily data from the grid

# Get weather data from the netCDF files for each year and each point of the
# grid that are within the study area.

years = np.arange(1961, 2017)
sta_in_zone = sta_gdf[sta_gdf['InZone'] == 1]
lat_indexes = list(sta_in_zone['lat_idx'])
lon_indexes = list(sta_in_zone['lon_idx'])
latitudes = list(sta_in_zone['latitude'])
longitudes = list(sta_in_zone['longitude'])

columns = pd.MultiIndex.from_tuples(
    list(zip(latitudes, longitudes)), names=['latitudes', 'longitudes'])
tasmax_dly = pd.DataFrame([], columns=columns)
tasmin_dly = pd.DataFrame([], columns=columns)
precip_dly = pd.DataFrame([], columns=columns)

for year in years:
    print("Fetching data for year %i... " % year, end='')
    filename = osp.join(infoclimat_dirname, "GCQ_v2_%i.nc" % year)
    netcdf_dset = netCDF4.Dataset(filename, 'r+')

    dates = pd.date_range(
        start=dt.datetime(year, 1, 1), end=dt.datetime(year, 12, 31))

    # Maximum dailty temperature.
    tasmax_dly = tasmax_dly.append(pd.DataFrame(
        np.array(netcdf_dset['tasmax'])[:, lat_indexes, lon_indexes],
        index=dates, columns=columns
        ))

    # Minimum daily temperature.
    tasmin_dly = tasmin_dly.append(pd.DataFrame(
        np.array(netcdf_dset['tasmin'])[:, lat_indexes, lon_indexes],
        index=dates, columns=columns
        ))

    # Total daily precipitation.
    precip_dly = precip_dly.append(pd.DataFrame(
        np.array(netcdf_dset['pr'])[:, lat_indexes, lon_indexes],
        index=dates, columns=columns
        ))

    netcdf_dset.close()
    print('done')

# %% Handle missing values

# print("Filling missing data... ", end='')

# # Fill missing daily precipitation with 0.
# precip_dly[precip_dly == -999] = 0

# # Fill missing daily air temperature with linear interpolation.
# tasmax_dly[tasmax_dly == -999] = np.nan
# tasmax_dly = tasmax_dly.interpolate(method='linear', axis=0)

# tasmin_dly[tasmin_dly == -999] = np.nan
# tasmin_dly = tasmin_dly.interpolate(method='linear', axis=0)

# print('done')

# %% Calculate monthly values

print("Calculating monthly values... ", end='')

tasmax_mly = tasmax_dly.groupby(
    [tasmax_dly.index.year, tasmax_dly.index.month]).mean()
tasmax_mly.index.rename(['Year', 'Month'], inplace=True)

tasmin_mly = tasmin_dly.groupby(
    [tasmin_dly.index.year, tasmin_dly.index.month]).mean()
tasmin_mly.index.rename(['Year', 'Month'], inplace=True)

precip_mly = precip_dly.groupby(
    [precip_dly.index.year, precip_dly.index.month]).sum()
precip_mly.index.rename(['Year', 'Month'], inplace=True)
precip_mly[precip_mly < 0] = -999

print('done')

# %% Save data to file

print("Saving data to file... ", end='')

tasmax_mly.to_csv(osp.join(
    project_dirname,
    'monthly_tasmax_pacc_({}-{}).csv'.format(
        np.min(tasmax_dly.index.year), np.max(tasmax_dly.index.year))
    ))
tasmin_mly.to_csv(osp.join(
    project_dirname,
    'monthly_tasmin_pacc_({}-{}).csv'.format(
        np.min(tasmin_dly.index.year), np.max(tasmin_dly.index.year))
    ))
precip_mly.to_csv(osp.join(
    project_dirname,
    'monthly_precip_pacc_({}-{})(missing).csv'.format(
        np.min(precip_dly.index.year), np.max(precip_dly.index.year))
    ))

print('done')
