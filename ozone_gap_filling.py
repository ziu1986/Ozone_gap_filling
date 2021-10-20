# Gap filling workflow
import os, glob, sys
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import matplotlib.pyplot as plt # Plotting
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
from scipy import fftpack
from scipy import stats
from bin.tools import *

def init():
    # Input directories
    src = os.environ['DATA']+'/astra/observations/ozone/'
    src_svanvik_OzoNorClim = src + '/Svanvik/NO0047R.*ozone*.xls'
    src_stations = ('Esrange', 'Jergul', 'Karasjok', 'Pallas', 'Svanvik')
    src_rra = os.environ['DATA']+'/nird/reanalysis/Copernicus/ensemble_ozone/SCA_ENSa.2018.O3.yearlyrea.nc'
    # Read data
    try:
        data
    except NameError:
        data = {}
        for station in src_stations:
            if station=='Barrow':
                data.update({station:load_data(src+station+'/*', type="Barrow")})
            else:
                data.update({station:load_data(src+station+'/*.nas')})
    # Concate Jergul and Karasjok data
    data.update({'jergkara':pd.concat((data['Jergul'], data['Karasjok']))})
    # Read and convert xls file data
    data_svanvik_OzoNorClim = []
    for file in sorted(glob.glob(src_svanvik_OzoNorClim)):
        tmp_data_svanvik = pd.read_excel(file, index_col=0, header=0)
        data_svanvik_OzoNorClim.append(tmp_data_svanvik['O3_mugm-3'].where(tmp_data_svanvik['O3_mugm-3']>=0.5).dropna()/2.)
    # Concat data Svanvik data
    data.update({'svanvik_OzoNorClim':pd.concat(data_svanvik_OzoNorClim)})
    # Load regional model reanalysis 2018 and set time axis
    try:
        data_rra = xr.open_dataset(src_rra)
        data_rra['time'] = pd.date_range("2018-01-01", periods=365*24, freq='H')
        data.update({'rra':data_rra})
    except NameError:
        print("Can't load regional data please check your source directory!")
    return(data)

def compute_time_lag(data):
    # Time lags -> fig6
    time_lag = range(-32,33)
    lag_jergkara_esrange = []
    lag_jergkara_pallas = []
    lag_svanvik_esrange = []
    lag_svanvik_pallas = []
    lag_svanvik_jergkara = []

    lag_label = ("jergkara_esrange","jergkara_pallas","svanvik_esrange","svanvik_pallas","svanvik_jergkara")
    for i in time_lag:
        lag_jergkara_esrange.append(time_lagged_corr(data['jergkara'], data['Esrange'], lag=i, pandas=True))
        lag_jergkara_pallas.append(time_lagged_corr(data['jergkara'], data['Pallas'], lag=i, pandas=True))
        lag_svanvik_esrange.append(time_lagged_corr(data['Svanvik'], data['Esrange'], lag=i, pandas=True))
        lag_svanvik_pallas.append(time_lagged_corr(data['Svanvik'], data['Pallas'], lag=i, pandas=True))
        lag_svanvik_jergkara.append(time_lagged_corr(data['Svanvik'], data['jergkara'], lag=i, pandas=True))

    # Print maximum in lag
    lag_max = {}
    print("Lag correlation")
    for i,lag in zip(lag_label,(lag_jergkara_esrange, lag_jergkara_pallas, lag_svanvik_esrange, lag_svanvik_pallas, lag_svanvik_jergkara)):
        lag_max.update({i:np.array(time_lag)[np.where(np.array(lag)==np.array(lag).max())[0]][0]})
        print("%s max at %d h" % (i, lag_max[i]))

    return(lag_max)

def main():
    print("Ozone gap filling")
    data = init()
    print("Loaded data: ", data.keys())
    lag_max = compute_time_lag(data)

if __name__ == "__main__":
    main()
