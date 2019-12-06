# -*- coding: utf-8 -*-
"""
Script to read and format NRCAN climate gridded data.

ftp://ftp.nrcan.gc.ca/pub/outgoing/canada_daily_grids/
"""

import csv
import numpy as np
import matplotlib.pyplot as plt
import re
import os.path as osp
import datetime as dt
from time import strftime


def read_daily_data(filename):
    """
    Fonction pour lire un feuillet de données météo de la grille NRCAN
    pour 1 journée.
    """
    with open(filename, 'r') as f:
        reader = list(csv.reader(f, delimiter=','))
        reader = reader[6:]

    data = [[]]
    for row in reader:
        values = re.findall('.{8}', row[0])
        data[-1].extend(values)
        if len(data[-1]) == 1068:
            data.append([])
    del data[-1]

    data = np.array(data).astype(float)
    data[data == -999] = np.nan

    return data


# %%
# Définir la position des centroïdes des mailles de la grille.

delta_deg = 300 / 3600
latitudes = np.array(
    list(reversed([41 + delta_deg/2 + delta_deg * i for i in range(510)])))
longitudes = np.array(
    [-141 + delta_deg/2 + delta_deg * i for i in range(1068)])


# %%
# Définir les indexes à utiliser pour générer une sous-grille.

idx_latitudes = np.where((latitudes >= 45) & (latitudes <= 47))[0]
idx_longitudes = np.where((longitudes >= -66) & (longitudes <= -64))[0]

latitudes = latitudes[idx_latitudes]
longitudes = longitudes[idx_longitudes]

# %%

data_stack = []
dirname = "D:/Data/grill_meteo_nrcan/canada_mintemperature"
for year in range(2009, 2011):
    for day in range(10):
        print('\rReading day {} of year {}   '.format(day + 1, year))
        filename = dirname + '/{}/min{}_{}.asc'.format(year, year, day + 1)
        if osp.exists(filename):
            data = read_daily_data(filename)
            data = data[idx_latitudes, :]
            data = data[:, idx_longitudes]
            data_stack.append(data)
print('')
data = np.stack(data_stack, axis=2)

_, _, nt = np.shape(data)
datetimes = [dt.datetime(2009, 1, 1) + dt.timedelta(days=1) * i for
             i in range(nt)]

# %%

npfname = ("D:/Data/grill_meteo_nrcan/data_tmin_sussex_2009-2017.npy")
content = {'data': data,
           'latitude': latitudes,
           'longitude': longitudes,
           'datetime': datetimes}
np.save(npfname, content, allow_pickle=True)

# %%

content = np.load(npfname, allow_pickle=True)

data = content.item()['data']
longitudes = content.item()['longitude']
latitudes = content.item()['latitude']
datetimes = content.item()['datetime']

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(datetimes, data[2, 2, :], '.-', lw=1,)
fig.autofmt_xdate()
fig.tight_layout()

fig2, ax2 = plt.subplots()
ax2.imshow(data[:, :, 5])
fig.tight_layout()
plt.show()

# %%

# Production des fichiers d'entrée pour PyHelp.

# Create an array of datestring and lat/lon
datestrings = [dt.strftime("%d/%m/%Y") for dt in datetimes]

ny, nx, nt = np.shape(data)
idx_row = np.repeat(np.arange(ny), nx)
idx_col = np.tile(np.arange(nx), ny)

# Reshape data into a 2d matrix.
data_2d = data[idx_row, idx_col, :].transpose()
lat_dd = latitudes[idx_row]
lon_dd = longitudes[idx_col]

# Remove cells with nan.
non_nan_idx = np.where(~np.isnan(data_2d[0, :]))[0]
data_2d = data_2d[:, non_nan_idx]
lat_dd = lat_dd[non_nan_idx]
lon_dd = lon_dd[non_nan_idx]

# Save the data to a csv file.
data_2d = data_2d.tolist()
lat_dd = lat_dd.tolist()
lon_dd = lon_dd.tolist()
varname = 'Minimum daily air temperature in \u00B0C'
fheader = [
    [varname],
    ['', ''],
    ['Created on ' + strftime("%d/%m/%Y")],
    ['Created from NRCAN grid'],
    ['', ''],
    ['Latitude (dd)'] + lat_dd,
    ['Longitude (dd)'] + lon_dd,
    ['', '']]

fdata = [[datestrings[i]] + data_2d[i] for i in range(nt)]
fcontent = fheader + fdata

fname = 'D:/Data/grill_meteo_nrcan/min_airtemp_input_data.csv'
with open(fname, 'w', encoding='utf8') as csvfile:
    writer = csv.writer(csvfile, delimiter=',', lineterminator='\n')
    writer.writerows(fcontent)
