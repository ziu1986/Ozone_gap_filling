# Tools for ozone_gap_filling
def datetime_from_time(start_date, time):
    '''
    Returns np.array of datetime objects from list of EMAC time (days since yyyy-mm-dd 00:00).
    Doesn't work for monthly mean.
    '''
    import numpy as np
    import datetime as dt

    # Getting hour from time
    delta_time = int(round((time[2]-time[1])*24))
    hour = (time-time.astype(int))*24 # np.array of floats
    day = time.astype(int) # np.array
    rHound = [] # round to full hour
    
    for iHour in hour:
        rHound.append(int(round(iHour)))
    minute = int(round((rHound[0]-hour[0])*60))
     # Generate the date vector
    dtime = np.array([ dt.timedelta(days=day[i], hours=hour[i], minutes=minute, seconds=0) for i in range(len(time)) ])
    x_time = start_date + dtime
    return x_time, delta_time

def read_station_data_noaa(infile, **kargs):
    '''
    Import of NOAA station data and return as pandas Series.
    **kargs utc - timeshift wrt UTC
            start_data - line number in which data block starts
            column - column index in which data is found
            station - standard, southpole, other format
    '''
    import codecs
    import datetime as dt
    import pandas as pd
    # Key words
    start_data = kargs.pop('start_data', 0)
    utc = kargs.pop('utc', 0)
    column = kargs.pop('column', 2)
    station = kargs.pop('station', 'standard')
    # Reading data
    lines = codecs.open(infile, "r", 'utf-8').read().splitlines()
    # Compute the timeshift
    timeshift = dt.timedelta(hours=utc)
    # Data structure
    time = []
    ozone_data_raw = []
    for each in lines[start_data:]:
        # Split the lines
        col = each.split()
        # Seperate ozone data
        ozone_data_raw.append(float(col[column]))
        # Generate timestamp
        if (station=='southpole') | (station=='other'):
            if  int(col[4]) == 24:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], '00'), '%Y-%m-%d %H')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s-%s-%s %s" % (col[1], col[2], col[3], col[4]), '%Y-%m-%d %H')
        else:
            if int(col[1][:2]) == 24:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" 00:00"), '%Y-%m-%d %H:%M')+dt.timedelta(days=1)
            else:
                temp_time = dt.datetime.strptime("%s" % (col[0]+" "+col[1]), '%Y-%m-%d %H:%M')
        time.append(temp_time-timeshift)
    # Mask the data (negative ozone values)
    ozone_data = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(ozone_data_raw))
    x_time = np.ma.masked_where(np.array(ozone_data_raw)<0, np.array(time))
    # Output pandas series
    pd_series = pd.Series(ozone_data, index=x_time)
    return pd_series


# Get the read-in routine for EBAS data
def read_station_data_ebas(infile, **kargs):
    '''
    Conversion of NASA aimes files from nilo.no data sets.
    '''
    import numpy as np
    import datetime as dt
    import pandas as pd
    import nappy as nap # Read and write NASA data files

    # Open file (read header only)
    nas_file = nap.openNAFile(infile)
    # Read actual data
    nas_file.readData()
    # Access the ozone data
    data_raw = np.array(nas_file.getNADict()['V'])
    # Close the nas-file
    nas_file.close()
    # Keywords
    conversion = kargs.pop('conversion', True)
    tracer = kargs.pop('tracer', ('O3',))
    verbose = kargs.pop('v', False)
    if (data_raw.shape[0]-1)/2 < len(tracer):
        for each in tracer:
            if nas_file.getNADict()['NCOM'][-1].find(each) > 0:
                return(read_station_data(infile, tracer=(each,)))
    # Conversion = 1/air_dens(25degC)/mass_fraction
    # 1/2 for O3 (0.5)
    # ca. 0.38 for SO2 
    # ca. 0.81 for NO
    M_air = 28.949                        # [g/mol] 
    air_dens_25 = 1.1839                  # [kg/m3]
    mass_fraction = {'O3':3*15.9994/M_air,
                     'SO2':32.065/M_air,# ug(S)/m3
                     'SO4':32.065/M_air,
                     'NO':14.0067/M_air,# ug(N)/m3
                     'NO2':14.0067/M_air}
       
    nas_date = nas_file.getNADict()['DATE']
    start_date = dt.datetime.strptime("%s-0%s-0%s 00:00:00" % (nas_date[0], nas_date[1], nas_date[2]), '%Y-%m-%d %H:%M:%S')
    # Filter the data (FLAG: 0 - valid, >0 - invalid)
    # Add day fraction at stop time to start date 
    x_time_station = np.ma.masked_where(data_raw[2]>0, 
                                        datetime_from_time(start_date, data_raw[0])[0])
    data_dic = {}
    if verbose:
        for each in tracer:
            print(each, data_raw.shape,
                  (nas_file.getNADict()['NCOM'][-1]).split().index(each),
                  (nas_file.getNADict()['NCOM'][-1]).split()[(nas_file.getNADict()['NCOM'][-1]).split().index(each)])
    #i = 1
    for each in tracer:
        # Find column in which tracer appears
        # Split table header and look for it
        # Table header has start/end time while
        # nas_file.getNADict()['V'] only gives one time field
        i = (nas_file.getNADict()['NCOM'][-1]).split().index(each)-1
        if conversion:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i]*1/air_dens_25/mass_fraction[each])
        else:
            data_dic[each] = np.ma.masked_where(data_raw[i+1]>0, data_raw[i])
        #i += 1
    pd_frame = pd.DataFrame(data_dic, index=x_time_station)
    return (pd_frame)

def load_data(src, **karg):
    '''
    Load ozone station data and round to full hours.
    Parameters
    ----------
    src : string
        Path to the files that shall be concatenated.
    Keyword arguments
    -----------------
    species : string
        Species that shall be extracted from file.
        Standard: "O3"
    type : string
        Type of source. Choices: ebas/barrow
    Returns
    -------
    pandas Timeseries.
    '''
    import os, glob, sys
    import pandas as pd
    
    species = karg.pop("species", "O3")
    src_type = karg.pop("type", "ebas")
    data = []
    for file in sorted(glob.glob(src)):
        print("Reading file %s" % (file))
        if src_type=="ebas":
            tmp = read_station_data_ebas(file)
            data.append(tmp[species]) 
        elif ((src_type=="Barrow") | (src_type=="barrow")) :
            if int(file[-4:]) < 2003:
                tmp = (read_station_data_noaa(file, utc=-9, start_data=28))
            elif int(file[-4:]) < 2012:
                tmp = (read_station_data_noaa(file, utc=-9, station='other', column=5))
            data.append(tmp)   
    # Concatenate the lists
    print('Concatenating data...')
    data = pd.concat(data)
    # Round to full hours
    data.index = data.index.round("h")
    return(data)

def time_lagged_corr(test_data, truth, **kargs):
    '''
    Compute time lagged correlation and returns correlation cooeficient.
    Takes test_data, truth as arguments. test_data is the one that gets shifted.
    kwargs:
       lag (0) - the time-lag in hours (can be positive or negative)
       v (False) - be a bit more verbose
       pandas (False) - use pandas dataframe object as input instead of numpy arrays
    '''
    import numpy as np
    import pandas as pd

    lag = kargs.pop('lag',0)
    verbose = kargs.pop('v', False)
    pandas = kargs.pop('pandas', False)
    if pandas:
        corr_coef = truth.corr(test_data.shift(lag))
        if verbose:
            corr_sign = corr_coef*np.sqrt((len(test_data)-np.fabs(lag)-2)/(1-corr_coef**2))
            print "%d %d %1.2f -> %2.2f" % (lag, len(test_data), corr_coef, corr_sign)
        return(corr_coef)
    else:
        if lag >= 0:
            corr_coef = np.ma.corrcoef(np.roll(test_data, lag)[lag:], truth[lag:])
        else:
            corr_coef = np.ma.corrcoef(np.roll(test_data, lag)[:lag], truth[:lag])
        corr_sign = corr_coef[0,1]*np.sqrt((len(test_data)-np.fabs(lag)-2)/(1-corr_coef[0,1]**2))
        if verbose:
            print "%d %d %1.2f -> %2.2f" % (lag, len(test_data), corr_coef[0,1], corr_sign)
        return corr_coef


def main():
    print("Tools for ozone_gap_filling")

if __name__ == "__main__":
    main()