import os, glob, sys, io
import numpy as np
import pandas as pd             # Timeseries data
import datetime as dt           # Time manipulation
import yaml
from matplotlib.dates import date2num # Convert dates to matplotlib axis coords
from matplotlib import dates
from scipy import fftpack
from scipy import stats
from bin.tools import *

def init():
    # Read configuration
    with open(r'config.yml') as file:
        config_list = yaml.load(file, Loader=yaml.FullLoader)
    src = config_list['sources']['ebas_ozone']
    src_svanvik_OzoNorClim = config_list['sources']['svanvik_ozone']
    src_rra = config_list['sources']['regional_ozone']
    station_list = config_list['station_list']
    file.close()
    # Read data
    try:    
        data = {}
        for station in station_list:
            if station=='Barrow':
                data.update({station:load_data(src+station+'/*', type="Barrow")})
            else:
                data.update({station:load_data(src+station+'/*.nas')})
    except NameError:
        sys.exit("Can't load ozone station data please check your source directory!")
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
        print("Warning: Can't load regional data please check your source directory!")
    return(data)

def extract_station_data(data, station_list):
    from bin.station_info import station_location
    local_rra = {}
    for each in station_list:
        local_rra.update({each:data['rra'].sel(lat=station_location[each].lat, lon=station_location[each].lon, method='nearest', time='2018-07')['O3']*0.5})
    return(local_rra)

def compute_time_lag(data):
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

def compute_clim(data):
    doys = np.arange(1,367)
    # Climatology from Esrange, Pallas, Jergul/Karasjok data
    climatology = pd.concat((data['Esrange'][:'2012'], data['Pallas'][:'2012'], data['jergkara'][:'2012']))

    # Daily mean climatology from Esrange, Pallas, Jergul/Karasjok data
    yozone, yerr, yerr_mean = compute_climatology(climatology)
    yozone_max, yerr_max, yerr_mean_max = compute_climatology(climatology, mode='max')
    yozone_min, yerr_min, yerr_mean_min = compute_climatology(climatology, mode='min')

    # Svanvik climatology
    yozone_svanvik, yerr_svanvik, yerr_mean_svanvik = compute_climatology(data['Svanvik'])
    yozone_max_svanvik, yerr_max_svanvik, yerr_mean_max_svanvik = compute_climatology(data['Svanvik'], mode='max')
    yozone_min_svanvik, yerr_min_svanvik, yerr_mean_min_svanvik = compute_climatology(data['Svanvik'], mode='min')

    # Hourly climatology
    clim_hourly, clim_hourly_err, clim_hourly_err_mean = compute_climatology(climatology, mode='hourly')
    clim_hourly_svanvik, clim_hourly_err_svanvik, clim_hourly_err_mean_svanvik = compute_climatology(data['Svanvik'], mode='hourly')

    # Compute spline fits
    from scipy.interpolate import UnivariateSpline
    # Fennoscandic climatology
    w = 1/yerr_mean
    fitSpl_dmean = UnivariateSpline(doys, climatology.groupby(climatology.index.dayofyear).apply(np.nanmean), w=w)
    dmax = climatology.resample('1d').apply(np.nanmax)
    fitSpl_dmax = UnivariateSpline(doys, dmax.groupby(dmax.index.dayofyear).apply(np.nanmean))
    # Svanvik
    w_svanvik = 1/yerr_mean_svanvik
    fitSpl_dmean_svanvik = UnivariateSpline(doys, data['Svanvik'].groupby(data['Svanvik'].index.dayofyear).apply(np.nanmean), w=w_svanvik)
    dmax_svanvik = data['Svanvik'].resample('1d').apply(np.nanmax)
    fitSpl_dmax_svanvik = UnivariateSpline(doys, dmax_svanvik.groupby(dmax_svanvik.index.dayofyear).apply(np.nanmean))
    
    # Pickle splines for comparison with other data
    import pickle
    with open('obs_climatologies.pkl','wb') as output:
        pickle.dump(fitSpl_dmean, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(fitSpl_dmean_svanvik, output, pickle.HIGHEST_PROTOCOL)
       
        pickle.dump(yerr_mean, output, pickle.HIGHEST_PROTOCOL)
        pickle.dump(yerr_mean_svanvik, output, pickle.HIGHEST_PROTOCOL)

    return({'clim':clim_hourly, 'clim_err':clim_hourly_err, 'clim_err_mean':clim_hourly_err_mean}, 
            {'clim':clim_hourly_svanvik, 'clim_err':clim_hourly_err_svanvik, 'clim_err_mean':clim_hourly_err_mean_svanvik})

def sample_climatology(clim, clim_svanvik):
    # Sample from houerly climatology
    sample_clim_svanvik = pd.DataFrame(pd.concat((clim_svanvik.iloc[:(31+28)*24],clim_svanvik.iloc[(31+29)*24:])).values, index=pd.date_range("2018-01-01 0:0", "2018-12-31 23:0", freq='H'))
    sample_clim = pd.DataFrame(pd.concat((clim.iloc[:(31+28)*24],clim.iloc[(31+29)*24:])).values, index=pd.date_range("2018-01-01 0:0", "2018-12-31 23:0", freq='H'))
    return(sample_clim, sample_clim_svanvik)

def compute_reconstruction(data, sample_clim, sample_clim_svanvik, lag_max, bias_corr = 1.2):
    # Bias correction for historical climatology to present day
    # Time lag correction (same for Esrange and Pallas)
    time_lag_corr = lag_max['svanvik_esrange']

    # Scaling factor
    scaling = sample_clim_svanvik/sample_clim.shift(-time_lag_corr)
    anomaly_pallas = data['Pallas']['07-2018']-sample_clim['07-2018'][0]
    anomaly_esrange = data['Esrange']['07-2018']-sample_clim['07-2018'][0]
    anomaly_svanvik = data['svanvik_OzoNorClim']['2018-07']-sample_clim_svanvik[0][data['svanvik_OzoNorClim']['2018-07'].index]-bias_corr

    reco_anomaly_svanvik = anomaly_pallas.shift(-time_lag_corr)*scaling['07-2018'][0]
    reco_svanvik = reco_anomaly_svanvik+sample_clim_svanvik['2018-07'][0]+bias_corr

    anomalies = {'Pallas': anomaly_pallas, 'Esrange':anomaly_esrange, 'Svanvik':anomaly_svanvik, 'Svanvik_reco':reco_anomaly_svanvik}

    return(anomalies, reco_svanvik, time_lag_corr)


def plot_reco(data, sample_clim, sample_clim_svanvik, anomalies, reco_svanvik, time_lag_corr):
    fig = plt.figure(1, figsize=(10,12))
    fig.canvas.set_window_title("ozone_reconstruction_2018_07")
    ax11 = plt.subplot(311)
    ax11.set_title('(a)')

    data['Esrange']['2018-07'].plot(ax=ax11, ls='None', marker='o', fillstyle='none', color='blue', label="Esrange")
    data['Pallas']['2018-07'].plot(ax=ax11, ls='None', marker='^', fillstyle='none', color='black', label="Pallas")
    data['svanvik_OzoNorClim']['2018-07'].plot(ax=ax11, color='blueviolet', ls='None', marker='d', label='Svanvik')
    sample_clim['07-2018'][0].plot(ax=ax11, color='red', label="Hourly clim.")
    sample_clim.shift(-time_lag_corr)['2018-07'][0].plot(ax=ax11, color='red', ls='--', label="Hourly clim. + time lag corr.")
    sample_clim_svanvik['07-2018'][0].plot(ax=ax11, color='red', ls='-.', label="Hourly clim. Svanvik")

    ax11.set_ylabel("$[O_3] (ppb)$")
    ax11.set_ylim(0,75)
    ax11.set_xticklabels("")
    ax11.set_xlabel('')
    ax11.legend(ncol=2)
    
    ax12 = plt.subplot(312)
    ax12.set_title('(b)')
    anomalies['Pallas'].plot(ax=ax12, ls='None', marker='^', fillstyle='none', color='black', label="Pallas")
    anomalies['Esrange'].plot(ax=ax12, ls='None', marker='o', fillstyle='none', color='blue', label="Esrange")
    anomalies['Svanvik'].plot(ax=ax12, ls='None', color='blueviolet', label='Svanvik', marker='d')
    anomalies['Svanvik_reco'].plot(ax=ax12, color='magenta', label='Reco. Svanvik')

    ax12.set_ylabel("$\Delta [O_3]$ (ppb)")
    ax12.set_ylim(-30, 30)
    ax12.set_xticklabels("")
    ax12.set_xlabel('')
    ax12.legend(ncol=2)

    ax13 = plt.subplot(313)
    ax13.set_title('(c)')

    reco_svanvik.plot(ax=ax13, ls='-', color='magenta', marker='None', label='Reco. Svanvik')
    data['svanvik_OzoNorClim']['2018-07'].plot(ax=ax13, color='blueviolet', fillstyle='none', ls='None', marker='d', label='Svanvik')
    try :
        data['svanvik_rra'].to_pandas().plot(ax=ax13, color='grey', fillstyle='none', ls=':', linewidth=2.5, label='CAMSRAQ')
    except KeyError:
        print("Warning: No regional data loaded!")

    ax13.set_ylabel("$[O_3] (ppb)$")
    ax13.set_ylim(0,75)
    ax13.set_xlabel('Time (days)')
    ax13.legend(ncol=3)

    return(fig)

def main():
    print("Tool kit for ozone gap filling")

if __name__ == "__main__":
    main()