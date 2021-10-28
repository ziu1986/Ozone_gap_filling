# Gap filling workflow
import matplotlib.pyplot as plt # Plotting
from bin.ogf_func import *

def main():
    print("Ozone gap filling")
    data = init()
    print("Loaded data: ", data.keys())
    lag_max = compute_time_lag(data)

    clim, clim_svanvik = compute_clim(data)
    print(clim.keys(), clim_svanvik.keys())

    sample_clim, sample_clim_svanvik = sample_climatology(clim['clim'], clim_svanvik['clim'])

    anomalies, reco_svanvik, time_lag_corr = compute_reconstruction(data, sample_clim, sample_clim_svanvik, lag_max)

    plot_reco(data, sample_clim, sample_clim_svanvik, anomalies, reco_svanvik, time_lag_corr)
    
    plot_all()

if __name__ == "__main__":
    main()
