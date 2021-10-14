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

def init():
    # Input directories
    src = os.environ['DATA']+'/astra/observations/ozone/'
    src_svanvik_OzoNorClim = src + '/Svanvik/NO0047R.*ozone*.xls'
    src_stations = ('Esrange', 'Janiskoski', 'Jergul', 'Karasjok', 'Pallas', 'Svanvik')
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
    except 
        print("Can'r load regional data please check your source directory!")
    return(data)

def main():
    print("Ozone gap filling")
    data = init()

if __name__ == "__main__":
    main()
