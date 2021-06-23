import datetime
import numpy as np
import pandas as pd
from matplotlib.dates import (date2num, DateFormatter)
from scipy import interpolate
from scipy.signal import lombscargle
from scipy.stats import pearsonr
import pickle
import matplotlib.pyplot as plt
from scipy.optimize import least_squares
from pathlib import Path
from os import listdir

# all code written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost
# run using python version 3.6

def glonasswlen(prn, signal):
    channel = [1, -4, 5, 6, 1, -4, 5, 6, -2, -7, 0, -1, -2, -7, 0, -1, 4, -3, 3, 2, 4, -3, 3, 2]
    try:
        channel_t = np.array([channel[i-1] for i in prn], dtype=float)
    except TypeError:
        channel_t = channel[prn-1]
    if signal == 'L1':
        lcar = 299792458 / (1602e06 + channel_t * 0.5625e06)
    elif signal == 'L2':
        lcar = 299792458. / (1246e06 + channel_t * 0.4375e06)
    else:
        print('signal not recognised')
        lcar = np.nan
    return lcar


def datetime2gps(dt):
    """
    converts a date to gps week and seconds of week
    not sure if it is better to use the matplotlib date as input?
    :param dt: in datetime format
    :return: gpsweek, sow (seconds of week) integers
    """
    dtgps = datetime.datetime(1980, 1, 6)
    deldtt = dt-dtgps
    gpsweek = int(deldtt.days/7)
    gpsweek_days = datetime.timedelta(days=gpsweek*7)
    sow = (dt - (dtgps+gpsweek_days))
    sow = sow.days*86400 + sow.seconds
    return gpsweek, sow


def ecef2lla(xyz):
    """
    convert earth centered earth fixed xyz position to lon, lat, altitude
    input should be in meters
    lon, lat in radians
    altitude in meters above WGS84 ellipsoid
    code translated from matlab function of same name written by Michael Kedler
    """
    # WGS84 constants
    a = 6378137
    e = 8.1819190842622e-2
    # calculations:
    b = np.sqrt(a**2 * (1 - e**2))
    ep = np.sqrt((a**2 - b**2) / b**2)
    p = np.sqrt(xyz[0]**2 + xyz[1]**2)
    th = np.arctan2(a*xyz[2], b*p)
    lon = np.arctan2(xyz[1], xyz[0])
    lat = np.arctan2((xyz[2] + ep**2*b*np.sin(th)**3), (p - e**2*a*np.cos(th)**3))
    N = a / np.sqrt(1 - e**2*np.sin(lat)**2)
    alt = p / np.cos(lat) - N
    # return lon in range[0, 2 * pi)
    lon = np.mod(lon, 2*np.pi)
    # correct for numerical instability in altitude near exact poles:
    # (after this correction, error is about 2 millimeters, which is about
    # the same as the numerical precision of the overall function)
    if np.abs(xyz[0]) < 1 and np.abs(xyz[1]) < 1:
        alt = np.abs(xyz[2]) - b
    return lat, lon, alt


def lla2ecef(lat, lon, alt):
    """
    convert lat, lon, altitude to earth centered earth fixed xyz position
    input in radians and meters above WGS84 ellipsoid
    code translated from matlab function of same name written by Michael Kedler
    """
    # WGS84 ellipsoid constants:
    a = 6378137
    e = 8.1819190842622e-2
    # intermediate calculation (prime vertical radius of curvature)
    N = a / np.sqrt(1 - e ** 2. * np.sin(lat) ** 2)
    # results:
    x = (N + alt) * np.cos(lat) * np.cos(lon)
    y = (N + alt) * np.cos(lat) * np.sin(lon)
    z = ((1 - e ** 2) * N + alt) * np.sin(lat)
    xyz = [x, y, z]
    return xyz


def gnss2azelv(xyz_ant, xyz_sat):
    """
    Converts position of antenna and position of satellite to
    azimuth and elevation angle of satellite relative to antenna
    :param xyz_sat: in meters
    :param xyz_ant: in meters
    :return: azi, elv in degrees
    """
    if len(xyz_sat) == 3:
        dx = xyz_sat[0] - xyz_ant[0]
        dy = xyz_sat[1] - xyz_ant[1]
        dz = xyz_sat[2] - xyz_ant[2]
    else:
        dx = xyz_sat[:, 0] - xyz_ant[0]
        dy = xyz_sat[:, 1] - xyz_ant[1]
        dz = xyz_sat[:, 2] - xyz_ant[2]
    lat, lon, _ = ecef2lla(xyz_ant)
    trans_matrix = np.zeros((3, 3))
    trans_matrix[0, 0] = -np.sin(lon)
    trans_matrix[1, 0] = -np.sin(lat) * np.cos(lon)
    trans_matrix[2, 0] = np.cos(lat) * np.cos(lon)
    trans_matrix[0, 1] = np.cos(lon)
    trans_matrix[1, 1] = -np.sin(lat) * np.sin(lon)
    trans_matrix[2, 1] = np.cos(lat) * np.sin(lon)
    trans_matrix[0, 2] = 0
    trans_matrix[1, 2] = np.cos(lat)
    trans_matrix[2, 2] = np.sin(lat)
    obsvec = np.vstack((dx, dy, dz))
    rss = np.matmul(trans_matrix, obsvec)
    azi = np.arctan2(rss[0], rss[1]) / np.pi * 180
    if any(azit < 0 for azit in azi):
        neginds = azi < 0
        azi[neginds] = 360 + azi[neginds]
    elv = np.arcsin(rss[2] / np.sqrt(np.sum(rss ** 2, axis=0))) / np.pi * 180
    return azi, elv


def readsp3file(sp3str):
    """
    Reads an sp3 orbit file and returns satellite identities and xyz coordinates
    input
    'sp3str': string to location of .sp3 orbit data
    outputs
    'sp3_df': pandas data frame object with 5 columns
    'DateTime' date and time of data in datetime.datetime format
    'sat_prn' sat constellation and number, e.g., 'G01' for GPS sat #1 ('R' for GLONASS, 'E' for Galileo)
    'X', 'Y' and 'Z' are ECEF coordinates (in m)
    """
    startid = 1
    counter = -1
    arrsize = (int(1440/5)+1)*92  # 5 mins plus one for midnight the next day, 92 satellites
    datetime_arr = np.empty(arrsize, dtype=object)
    satprn_arr = np.empty(arrsize, dtype=object)
    X_arr = np.empty(arrsize, dtype=float)
    Y_arr = np.empty(arrsize, dtype=float)
    Z_arr = np.empty(arrsize, dtype=float)
    with open(sp3str, 'r') as f:
        for line in f:
            if line[0:3] == 'EOF':
                break
            if line[0] == '*':
                curdt = datetime.datetime(int(line[3: 7]), int(line[8: 10]), int(line[11: 13]),
                                          int(line[14: 16]), int(line[17: 19]), int(line[20: 22]))
                startid = 0
            elif startid == 0 and line[0] != '*':
                # can get rid of this condition if you care about other satellites
                if line[1:2] in ['G', 'R', 'E']:
                    counter = counter + 1
                    satprn_arr[counter] = line[1:4]
                else:
                    continue
                datetime_arr[counter] = curdt
                X_arr[counter] = float(line[4: 18]) * 1000
                Y_arr[counter] = float(line[18: 32]) * 1000
                Z_arr[counter] = float(line[32: 46]) * 1000
    if counter+1 < arrsize:
        datetime_arr = datetime_arr[:counter]
        satprn_arr = satprn_arr[:counter]
        X_arr = X_arr[:counter]
        Y_arr = Y_arr[:counter]
        Z_arr = Z_arr[:counter]
    sp3_data = {'DateTime': datetime_arr, 'sat_prn': satprn_arr, 'X': X_arr, 'Y': Y_arr, 'Z': Z_arr}
    sp3_df = pd.DataFrame(sp3_data)
    return sp3_df


def nmea2snr(nmeastr, snrdir, sp3dir=False, **kwargs):
    """
    This function reads NMEA 0183 format GPS data (e.g., recorded by low-cost GNSS hardware) and
    converts it into organised files for GNSS-R analysis
    :param nmeastr: path to nmea file
    :param snrdir: path to directory where organised SNR data is saved as 'pickle' files
    :param sp3dir: path to directory with sp3 orbit data,
    if sp3dir=False then saves the azimuth and elevation angle data from NMEA file
    to load use the following code:
    f = open(snrstr, 'rb')
    snr_df = pickle.load(f)  # SNR data
    fix_df = pickle.load(f)  # fix data (lat, lon, height)
    f.close()
    OUTPUT DATA FORMAT
    'snr_df' is in pandas.DataFrame format with the following columns
    'datenum': time of observation in maplotlib dates format
    'sat_prn': satellite constellation and number in the same format as in '.SP3' files (e.g., G01 = GPS sat 1)
    'elevation': elevation angle of satellite in degrees
    'azimuth': azimuth angle of satellite in degrees
    'snr': SNR in units of dB-Hz
    'seconds': seconds of day of observation
    'fix_df' is in pandas.DataFrame format with the following columns
    'seconds': seconds of day
    'latitude': in degrees
    'longitude': in degrees
    'height': in m relative to WGS84 elipsoid
    :**kwargs:
    elvlims: elevation angle limits (default [0, 90])
    azilims: azimuth angle limits (default = [0, 360])
    tempres: adjust temporal resolution in seconds (default is resolution of nmea file)
    xyz_ant: earth centered earth fixed xyz coordinates of antenna
    (if xyz_ant not given then it will be taken as average fix data)
    """
    elvlims = kwargs.pop('elvlims', [0, 90])
    azilims = kwargs.pop('azilims', [0, 360])
    donefile = 0
    cnt = 0
    days_all = np.empty(0, dtype=object)
    while donefile == 0:
        with open(nmeastr, 'r') as f:
            line = f.readline()
            fixdata = np.empty((86400, 4), dtype=float)
            fixdata[:] = np.nan
            nmeadata = {}
            if cnt != 0:
                cnt2 = 0
                while cnt2 < cnt:
                    cnt2 = cnt2 + 1
                    line = f.readline()
                    continue
            print('sorting through nmea obs ...')
            sect = np.nan
            while line:
                row = line.split(',')
                if '$' in line[1:-1]:
                    print('misplaced $ issue')
                    line = f.readline()
                    cnt = cnt + 1
                    continue
                if row[0][3:6] == 'RMC':
                    sect = int(float(row[1][0:2]) * 60 * 60 + float(row[1][2:4]) * 60 + float(row[1][4:]))
                    curdt = datetime.datetime(int(row[9][4:6])+2000,
                                              int(row[9][2:4]),
                                              int(row[9][0:2]))
                    if curdt not in days_all:
                        days_all = np.append(days_all, curdt)
                        curdt = curdt - datetime.timedelta(days=1)
                        if len(days_all) > 1:
                            print('got one day of data:')
                            print(curdt)
                            break
                elif row[0][3:6] == 'GGA':
                    sect = int(float(row[1][0:2]) * 60 * 60 + float(row[1][2:4]) * 60 + float(row[1][4:]))
                    latt = float(row[2][0:2]) + float(row[2][2:])/60
                    if row[3] == 'S':
                        latt = -latt
                    lont = float(row[4][0:3]) + float(row[4][3:])/60
                    if row[5] == 'W':
                        lont = -lont
                    hgtt = float(row[9])+float(row[11])
                    fixdata[sect, :] = [sect, latt, lont, hgtt]
                elif row[0][3:6] == 'GSV' and np.isnan(sect) == 0:
                    satconst_nmea = row[0][2]
                    if satconst_nmea == 'P':
                        satconst_sp3 = 'G'  # GPS
                    elif satconst_nmea == 'L':
                        satconst_sp3 = 'R'  # GLONASS
                    elif satconst_nmea == 'A':
                        satconst_sp3 = 'E'  # GALILEO
                    else:
                        satconst_sp3 = satconst_nmea
                    ind = 4  # starting from 4th entry
                    while ind < len(row) - 2:
                        lens = [len(row[ind]), len(row[ind + 1]), len(row[ind + 2]), len(row[ind + 3])]
                        if all(lent != 0 for lent in lens):
                            satoffset = 0
                            if satconst_sp3 == 'G' and int(row[ind]) > 32:
                                ind = ind + 4
                                continue
                            if satconst_sp3 == 'R':
                                satoffset = 64
                            satnostr = str(int(row[ind]) - satoffset)
                            if len(satnostr) == 1:
                                satnostr = '0' + satnostr
                            satidt = satconst_sp3 + satnostr
                            if satidt not in nmeadata:
                                nmeadata[satidt] = np.empty((86400, 4))
                                nmeadata[satidt][:] = np.nan
                            nmeadata[satidt][sect, :] = [sect, int(row[ind + 1]), int(row[ind + 2]), int(row[ind + 3])]
                            # elevation, azimuth, SNR
                        ind = ind + 4
                line = f.readline()
                cnt = cnt + 1
        if not line:
            donefile = 1
            print('end of file')
        if cnt == 0:
            print('empty file')
            return
        f.close()

        if 'xyz_ant' in kwargs:
            xyz_ant = kwargs.get('xyz_ant')
        else:
            # now calculate guess of xyz_ant
            meanlat = np.nanmedian(fixdata[:, 1])
            meanlon = np.nanmedian(fixdata[:, 2])
            meanhgt = np.nanmedian(fixdata[:, 3])
            # print mean values if you want
            # print('mean lat is ' + str(meanlat))
            # print('mean lon is ' + str(meanlon))
            # print('mean hgt is ' + str(meanhgt))
            # then use lla2ecef, not the end of the world if not too precise
            xyz_ant = lla2ecef(meanlat/180*np.pi, meanlon/180*np.pi, meanhgt)

        dt = 1  # maximum possible resolution
        if 'tempres' in kwargs:
            dt = kwargs.get('tempres')

        if sp3dir:
            # get date strings for finding the sp3 files
            doy_str = str(curdt.timetuple().tm_yday)
            year_str = str(curdt.year)
            # you could add more options for the orbit files if you want
            # since the naming conventions are different, give diff types of files unique directory names
            if sp3dir[-3:] == 'COD':
                sp3str = sp3dir + '/COD0MGXFIN_' + year_str + doy_str + '0000_01D_05M_ORB.SP3'
            elif sp3dir[-3:] == 'GFZ':
                sp3str = sp3dir + '/GFZ0MGXRAP_' + year_str + doy_str + '0000_01D_05M_ORB.SP3'
            else:
                print('you need to configure this code for different SP3 orbit files')
                return
            print('getting orbit info for ' + str(curdt.date()))
            sp3_df = readsp3file(sp3str)
        else:
            print('no sp3 data - using azimuth and elevation values directly from NMEA')

        # now doing putting into snr_data format
        snrdata = np.empty((0, 6))
        for satidt in nmeadata:
            nmeadatat = nmeadata[satidt]
            # first get rid of nans
            tfilter = np.isnan(nmeadatat[:, 0]) == 0
            nmeadatat = nmeadatat[tfilter, :]
            # now do the dt adjustment
            tfilter = np.mod(nmeadatat[:, 0], dt) == 0
            nmeadatat = nmeadatat[tfilter, :]
            # now make a datetime array
            datetime_t = curdt + nmeadatat[:, 0] * datetime.timedelta(seconds=1)
            # now make a satellite a
            tempsats = np.empty((len(datetime_t)), dtype=object)
            tempsats[:] = satidt
            tempdata = np.column_stack((date2num(datetime_t), tempsats, nmeadatat[:, 1], nmeadatat[:, 2],
                                        nmeadatat[:, 3], nmeadatat[:, 0]))
            if sp3dir and satidt in sp3_df['sat_prn'].values:
                tfilter = sp3_df['sat_prn'] == satidt
                tt_sp3 = sp3_df['DateTime'][tfilter].values
                tt_sp3 = [(tt - tt_sp3[0]).astype('timedelta64[s]').astype(int) for tt in tt_sp3]
                tt_sp3_new = np.linspace(tt_sp3[0], tt_sp3[-1], int((tt_sp3[-1] - tt_sp3[0]) / dt) + 1)
                t_sp3_new = np.array(tt_sp3_new, dtype=float)
                xt = sp3_df['X'][tfilter].values
                xtck = interpolate.splrep(tt_sp3, xt)
                x_new = interpolate.splev(tt_sp3_new, xtck)
                yt = sp3_df['Y'][tfilter].values
                ytck = interpolate.splrep(tt_sp3, yt)
                y_new = interpolate.splev(tt_sp3_new, ytck)
                zt = sp3_df['Z'][tfilter].values
                ztck = interpolate.splrep(tt_sp3, zt)
                z_new = interpolate.splev(tt_sp3_new, ztck)
                xyz_new = np.column_stack((x_new, y_new, z_new))
                azit, elvt = gnss2azelv(xyz_ant, xyz_new)
                # now find the overlapping dates
                _, ind_nmea, ind_sp3 = np.intersect1d(np.array(nmeadatat[:, 0], dtype=float), tt_sp3_new,
                                                      return_indices=True)
                tempdata = tempdata[ind_nmea]
                tempdata[:, 2] = elvt[ind_sp3]
                tempdata[:, 3] = azit[ind_sp3]
            elif sp3dir:
                print('missing orbit data for ' + satidt)
                tempdata = []
            # then collect all the data
            if len(tempdata) > 0:
                snrdata = np.vstack((snrdata, tempdata))
        tfilter = np.logical_and(snrdata[:, 2] >= elvlims[0], snrdata[:, 2] <= elvlims[1])
        snrdata = snrdata[tfilter]
        tfilter = np.logical_and(snrdata[:, 3] >= azilims[0], snrdata[:, 3] <= azilims[1])
        snrdata = snrdata[tfilter]
        snr_df = pd.DataFrame({'datenum': snrdata[:, 0], 'sat_prn': snrdata[:, 1], 'elevation': snrdata[:, 2],
                               'azimuth': snrdata[:, 3], 'snr': snrdata[:, 4], 'seconds': snrdata[:, 5]})
        nanfilter = np.isnan(fixdata[:, 0]) == 0
        fixdata = fixdata[nanfilter, :]
        fix_df = pd.DataFrame({'seconds': fixdata[:, 0], 'latitude': fixdata[:, 1], 'longitude': fixdata[:, 2],
                               'height': fixdata[:, 3]})
        # the 'seconds' is mainly for finding unique data more easily
        snr_df = snr_df.sort_values(by=['datenum'], ignore_index=True)
        # now need to pickle and save data
        Path(snrdir).mkdir(parents=True, exist_ok=True)
        # first check if data already exists:
        snrfilestr = snrdir + '/' + str(curdt.strftime("%y_%m_%d")) + '.pkl'
        try:
            f = open(snrfilestr, 'rb')
            snr_df_old = pickle.load(f)
            frames = [snr_df_old, snr_df]
            snr_df = pd.concat(frames, ignore_index=True)
            snr_df = snr_df.drop_duplicates(subset=['seconds', 'sat_prn'])
            snr_df = snr_df.sort_values(by=['datenum'], ignore_index=True)
            fix_df_old = pickle.load(f)
            frames = [fix_df_old, fix_df]
            fix_df = pd.concat(frames, ignore_index=True)
            fix_df = fix_df.drop_duplicates(subset=['seconds'])
            fix_df = fix_df.sort_values(by=['seconds'], ignore_index=True)
            f.close()
        except IOError:
            print('new file')
        f = open(snrfilestr, 'wb')
        pickle.dump(snr_df, f)
        pickle.dump(fix_df, f)
        f.close()
        print('dumped a pickle')
    return


def snr2arcs(snr_df, rhlims, signal='L1', arctlim=False, pktnlim=0, normalize=False, snrfigs=False, lspfigs=False,
             polydeg=2, **kwargs):
    """
    reads a formatted snr DataFrame (e.g., from nmea2snr) and organises into:
    reflector height estimates (+other stats) and detrended snr arcs for inverse analysis
    :param snr_df: pandas data frame of organised SNR data
    :param rhlims: upper and lower reflector height limits (in metres) for quality control
    :param signal: default 'L1' (C/A), can also use L2...if want to use L2C or L5, need to make some edits
    :param arctlim: split arcs into subarcs of time in units of seconds (default = False)
    :param pktnlim: QC control if you want to use a minimum peak-to-noise limit for periodograms (default = 0)
    :param normalize: if you want to normalize the arcs so that they have the same amplitude (default = False)
    :param snrfigs: if you want to produce some figures of SNR arcs (default = True)
    :param lspfigs: if you want to produce some figures of Lomb-Scargle Periodograms (default = False)
    :param polydeg: degree of polynomial to fit and subtract from SNR data
    :param kwargs: see below
    elvlims: elevation angle limits (e.g., [5, 30])
    azilims: azimuth angle limits (e.g., [90, 270])
    tempres: if want to use different temporal resolution to input data
    satconsts: default use all given, otherwise specify from ['G', 'R', 'E'] (gps / glonass / galileo)
    :return rh_df: dataframe of reflector height estimates and stats
    :return snr_dt: detrended SNR arcs e.g., for inverse analysis
    """

    # now need to convert SNR to linear scale (from dB-Hz)
    snr_df['snr'] = 10 ** (snr_df['snr'].values / 20)

    if 'tempres' in kwargs:
        tempres = kwargs.get('tempres')
        # now do the modulate
        tfilter = np.mod(snr_df['seconds'].values, tempres) == 0
        snr_df = snr_df[tfilter]
    else:
        # find the temporal resolution
        secst = np.unique(snr_df['seconds'])
        tempres = secst[1] - secst[0]
    # print('temporal resolution is ' + str(int(tempres)) + ' seconds')

    if 'satconsts' in kwargs:
        satconsts = kwargs.get('satconsts')
        # snr_df['sat_prn'].astype(str).str[0]
        allsats = ['G', 'R', 'E']
        for satc in allsats:
            if satc not in satconsts:
                tfilter = snr_df['sat_prn'].astype(str).str[0] != satc
                snr_df = snr_df[tfilter]

    # elvlims and azilims
    if 'elvlims' in kwargs:
        elvlims = kwargs.get('elvlims')
        tfilter = np.logical_and(snr_df['elevation'] > elvlims[0], snr_df['elevation'] < elvlims[1])
        snr_df = snr_df[tfilter]
    if 'azilims' in kwargs:
        azilims = kwargs.get('azilims')
        tfilter = np.logical_and(snr_df['azimuth'] > azilims[0], snr_df['azimuth'] < azilims[1])
        snr_df = snr_df[tfilter]

    rh_arr = np.empty((0, 14))
    snr_dt = np.empty((0, 4), dtype=object)
    for sat in np.unique(snr_df['sat_prn']):
        tfilter = snr_df['sat_prn'] == sat
        temp_df = snr_df[tfilter]
        if sat[0] == 'G' or sat[0] == 'E':
            if signal == 'L1':
                lcar = 299792458 / 1575.42e06
            elif signal == 'L2':
                lcar = 299792458 / 1227.60e06
            else:
                print('signal not recognised')
        elif sat[0] == 'R':
            lcar = glonasswlen(int(sat[1:3]), signal)
        maxf = 2 * (rhlims[0] + rhlims[1]) / lcar
        precisionf = 2 * 0.001 / lcar  # 1 mm
        f = np.linspace(precisionf, maxf, int(maxf / precisionf))
        # then look for points where either
        # 1 there is a gap in data bigger 2 or 3 times the temporal resolution
        # or where the elevation changes from ascending to descending
        temp_df = temp_df.sort_values(by='seconds', ignore_index=True)
        ddate = np.ediff1d(temp_df['datenum'])
        bkpt = np.where(ddate > 1/(24*12))[0]  # gaps bigger than 5 mins
        delv = np.ediff1d(temp_df['elevation'])
        bkpt = np.append(bkpt, np.where(np.diff(np.sign(delv)))[0])
        bkpt = np.unique(bkpt)
        bkpt = np.sort(bkpt)
        cnt = -1
        while cnt < len(bkpt):
            cnt = cnt + 1
            if cnt == 0:
                sind = 0
            else:
                sind = bkpt[cnt-1] + 1
            if cnt == len(bkpt):
                eind = len(delv)
            else:
                eind = bkpt[cnt] + 1
            if eind-sind < 21:
                continue
            elvt = np.array(temp_df['elevation'][sind:eind].values, dtype=float)
            sinelvt = np.sin(elvt / 180 * np.pi)
            snrt = np.array(temp_df['snr'][sind:eind].values, dtype=float)
            datet = np.array(temp_df['datenum'][sind:eind].values, dtype=float)
            azit = np.array(temp_df['azimuth'][sind:eind].values, dtype=float)
            z = np.polyfit(sinelvt, snrt, polydeg)
            p = np.poly1d(z)
            snrdt = snrt - p(sinelvt)
            # start doing subarcs here
            sindt = 0
            cnt2 = 0
            while sindt < len(sinelvt):
                if not arctlim:
                    eindt = len(sinelvt)
                elif datet[sindt] + (arctlim-tempres/2)/86400 < datet[-1]:
                    sdatet = datet[sindt]
                    edatet = sdatet + (arctlim-tempres/2)/86400
                    eindt = np.min(np.where(datet > edatet)) + 1  # remember the +1 for stupid indexing
                else:
                    break
                if eindt-sindt < 21:
                    break
                # add a length / time condition
                sinelvtt = sinelvt[sindt:eindt]
                snrdtt = snrdt[sindt:eindt]
                if normalize:
                    snrdtt = snrdtt * 100 / (np.max(np.abs(snrdtt)))
                datett = datet[sindt:eindt]
                elvtt = elvt[sindt:eindt]
                azitt = azit[sindt:eindt]
                sattt = np.empty((len(datett)), dtype=object)
                sattt[:] = sat
                temp_arr = np.column_stack([datett, sattt, sinelvtt, snrdtt])
                snr_dt = np.vstack((snr_dt, temp_arr))
                pgram = lombscargle(sinelvtt, snrdtt, f * 2 * np.pi, normalize=True)
                reflh = 0.5 * f * lcar
                tfilter = np.logical_and(reflh > rhlims[0], reflh < rhlims[1])
                pgram_sub = pgram[tfilter]
                reflh_sub = reflh[tfilter]
                maxind = np.argmax(pgram_sub)
                tfilter = np.logical_or(reflh < rhlims[0], reflh > rhlims[1])
                pgram_outl = pgram[tfilter]
                pktn = np.max(pgram_sub) / np.mean(pgram_outl)
                # if maxind == 0 or maxind == len(pgram_sub)-1:
                #     print('!!!!!!!!!!!!!!!!!!!!COUNT!!!!!!!!!!!!!!!!!!!!')
                if pktn > pktnlim and maxind != 0 and maxind != len(pgram_sub)-1:
                    cnt2 = cnt2 + 1
                    temp_arr = np.empty((1, 14), dtype=object)
                    temp_arr[0, 0] = np.mean(datett)  # time of arc
                    temp_arr[0, 1] = reflh_sub[maxind]
                    temp_arr[0, 2] = sat
                    dthdt = ((elvtt[-1] - elvt[0]) / 180 * np.pi) / (86400 * (datett[-1] - datett[0]))
                    temp_arr[0, 3] = np.tan(np.mean(elvtt) / 180 * np.pi) / dthdt
                    temp_arr[0, 4] = np.min(elvtt)
                    temp_arr[0, 5] = np.max(elvtt)
                    temp_arr[0, 6] = np.mean(azitt)
                    temp_arr[0, 7] = np.mean(p(sinelvtt))
                    temp_arr[0, 8] = np.max(pgram_sub)
                    temp_arr[0, 9] = np.var(snrdtt)
                    temp_arr[0, 10] = 86400 * (datett[-1] - datett[0])
                    temp_arr[0, 11] = pktn
                    temp_arr[0, 12] = len(elvtt)
                    temp_arr[0, 13] = len(f)
                    rh_arr = np.vstack((rh_arr, temp_arr))
                    if lspfigs:
                        fig, ax = plt.subplots(figsize=(6.5, 2.5))
                        ax.plot(reflh, pgram)
                        plt.savefig(sat + '_' + str(cnt) + str(cnt2) + '_LSP' + '.png')
                        plt.close(fig)
                    if snrfigs:
                        fig, ax = plt.subplots(figsize=(6.5, 2.5))
                        ax.plot(elvtt, snrdtt)
                        plt.xlim(elvlims[0], elvlims[1])
                        plt.ylim(-220, 220)
                        plt.savefig(sat + '_' + str(cnt) + str(cnt2) + '_SNR' + '.png')
                        plt.close(fig)
                sindt = eindt
    rh_df = pd.DataFrame({'datenum': rh_arr[:, 0], 'refl_hgt': rh_arr[:, 1], 'sat_prn': rh_arr[:, 2],
                          'tane_dedt': rh_arr[:, 3], 'min_e': rh_arr[:, 4], 'max_e': rh_arr[:, 5],
                          'mean_azi': rh_arr[:, 6], 'mean_power': rh_arr[:, 7], 'pgram_peak': rh_arr[:, 8],
                          'variance': rh_arr[:, 9], 'arc_length': rh_arr[:, 10], 'pktn_ratio': rh_arr[:, 11],
                          'numpts': rh_arr[:, 12], 'numfs': rh_arr[:, 13]})
    rh_df = rh_df.sort_values(by=['datenum'], ignore_index=True)
    snrdt_df = pd.DataFrame({'datenum': snr_dt[:, 0], 'sat_prn': snr_dt[:, 1],
                             'sin_e': snr_dt[:, 2], 'snr_dt': snr_dt[:, 3]})
    snrdt_df = snrdt_df.sort_values(by=['datenum'], ignore_index=True)
    return rh_df, snrdt_df


def residuals_bsp_spectral(sfacs, bspline_order, knots, rh_df):
    dt_even = 1 / (24*60)
    t_even = np.linspace(knots[0], knots[-1], int((knots[-1] - knots[0]) / dt_even))
    bspl_even = interpolate.splev(t_even, (knots, sfacs, bspline_order))
    dhdt_even = np.gradient(bspl_even, dt_even)
    f = interpolate.interp1d(t_even, dhdt_even)
    datenumt = np.array(rh_df['datenum'].values, dtype=float)
    dhdt = f(datenumt)
    dhdt = dhdt/86400
    tane_dedt = rh_df['tane_dedt'].values
    bspl_adj = interpolate.splev(datenumt, (knots, sfacs, bspline_order)) + dhdt * tane_dedt
    residual_spectral = np.array(bspl_adj - rh_df['refl_hgt'].values, dtype=float)
    return residual_spectral


def residuals_bsp_js(sfacs, bspline_order, knots, satconsts, signal, snrdt_df):
    if len(sfacs) - len(satconsts) * 2 == len(knots) - bspline_order - 1:
        # then no roughness
        hgt_sfacs = sfacs[: - len(satconsts) * 2]
        satparams = sfacs[- len(satconsts) * 2:]
        rough_in = False
    elif len(sfacs) - len(satconsts) * 2 - 1 == len(knots) - bspline_order - 1:
        # roughness
        hgt_sfacs = sfacs[: - len(satconsts) * 2 -1]
        satparams = sfacs[- len(satconsts) * 2 -1: -1]
        rough_in = sfacs[-1]
    else:
        print('issue with length of input parameter array')
        return
    tmpc = 0
    res = np.empty(0)
    for satc in satconsts:
        tmpc = tmpc + 1
        tfilter = snrdt_df['sat_prn'].astype(str).str[0] == satc
        tmp_df = snrdt_df[tfilter]
        if satc == 'G' or satc == 'E':
            if signal == 'L1':
                lcar = 299792458 / 1575.42e06
            elif signal == 'L2':
                lcar = 299792458 / 1227.60e06
            else:
                print('signal not recognised')
        elif satc == 'R':
            satnos = np.array(tmp_df['sat_prn'].astype('string').str[1:3], dtype=int)
            lcar = glonasswlen(satnos, signal)
        datenumt = np.array(tmp_df['datenum'].values, dtype=float)
        sin_et = np.array(tmp_df['sin_e'].values, dtype=float)
        snr_dtt = np.array(tmp_df['snr_dt'].values, dtype=float)
        h1 = interpolate.splev(datenumt, (knots, hgt_sfacs, bspline_order))
        modelsnr = satparams[int((tmpc - 1)*2)]  * np.sin(4 * np.pi * h1 * sin_et / lcar) + \
                   satparams[int((tmpc - 1)*2 + 1)] * np.cos(4 * np.pi * h1 * sin_et / lcar)
        if rough_in:
            lk = 2 * np.pi / lcar
            modelsnr = modelsnr * np.exp(-4 * lk ** 2 * rough_in * sin_et ** 2)
        res = np.append(res, modelsnr - snr_dtt)
    return res


def invsnr(sdatetime, edatetime, snrdir, invdir, kspac, tlen, rhlims, snrfit=True, signal='L1', largetides=True,
           arctlim=False, pktnlim=0, normalize=False, smoothqc=True, bspline_order=2, **kwargs):
    """

    :param invdir: string to directory where output files are saved as '.pkl' files, then plotted using 'invsnr_plot'
    :param sdatetime: datetime format of start date and time
    :param edatetime: datetime format of end date and time
    :param snrdir: string to directory where snr data is contained (i.e., from nmea2snr function output)
    :param kspac: knot spacing - see Strandberg et al. (2016) or other b-spline literature
    :param tlen: time window length - see Strandberg et al. (2016), tlen should be at least 3 x bigger than kspac.
    tlen should be a multiple of 3 hours - the algorithm moves forwards in steps of tlen/3
    :param rhlims: reflector height limits, e.g., [4, 6] (in meters)
    :param snrfit: set to False if you don't want to do the full inverse modelling
    :param signal: 'L1' or 'L2'
    :param largetides: set to False if you'd rather not guess the b-spline scaling factors prior to iversion
    :param arctlim: split satellite arcs into sub arcs of length arctlim (in seconds)
    :param pktnlim: QC condition, ratio of periodogram peak to noise must be greater than this value
    :param normalize: normalize the snr prior to inversion
    :param smoothqc: qc to filter out refl_hgt values > 3 std away from smoothed signal using mov_avg
    :param bspline_order: default = 2, not worth going higher and definitely don't want 1 (linear)
    :param kwargs:
    elvlims: elevation angle limits (e.g., [5, 30])
    azilims: azimuth angle limits (e.g., [90, 270])
    tempres: if want to use different temporal resolution to input data
    satconsts: default use all given, otherwise specify from ['G', 'R', 'E'] (gps / glonass / galileo)
    rough_in: guess of sea surface roughness to fit snr data, see Strandberg et al, (2016)
    :return:
    """

    tlen_td = datetime.timedelta(hours=int(tlen*24))
    tdatetime = sdatetime - 2*tlen_td/3
    while tdatetime < edatetime - 2*tlen_td/3:
        invout = {}
        tdatetime = tdatetime + tlen_td/3
        tdatetime_end = tdatetime + tlen_td
        tdatenum = date2num(tdatetime)
        tdatenum_end = date2num(tdatetime_end)
        curdt = tdatetime + tlen_td/3
        print(str(curdt.date()) + ' ' + str(curdt.time()))
        # first work out how many days of data need to load
        mlen = 1
        if tdatetime.day != tdatetime_end.day:
            # mlen = tdatetime_end.day - tdatetime.day + 1
            mlen = tdatenum_end - np.mod(tdatenum_end, 1) - (tdatenum - np.mod(tdatenum, 1)) + 1
            if tdatetime_end.hour == 0:
                mlen = mlen - 1
        mlen = int(mlen)
        snr_df = pd.DataFrame({})
        missingdata = False
        for m in range(1, int(mlen)+1):
            curdtt = tdatetime + datetime.timedelta(days=m-1)
            snrfilestr = snrdir + '/' + str(curdtt.strftime("%y_%m_%d")) + '.pkl'
            try:
                f = open(snrfilestr, 'rb')
                snr_dft = pickle.load(f)
                f.close()
            except IOError:
                print('missing data')
                missingdata = True
                break
            tfilter = np.logical_and(snr_dft['datenum'] >= tdatenum,
                                     snr_dft['datenum'] < tdatenum_end)
            snr_dft = snr_dft[tfilter]
            frames = [snr_df, snr_dft]
            snr_df = pd.concat(frames, ignore_index=True)
            if len(snr_df.index) == 0:
                print('no data')
                missingdata = True
                break
        if missingdata:
            continue
        rh_df, snrdt_df = snr2arcs(snr_df, rhlims, signal=signal, arctlim=arctlim, pktnlim=pktnlim, normalize=normalize,
                                   **kwargs)

        if 'satconsts' in kwargs:
            satconsts = kwargs.get('satconsts')
        else:
            satconsts = rh_df['sat_prn'].astype(str).str[0]
            satconsts = np.unique(satconsts)
        if np.min(rh_df['datenum'].values) > date2num(tdatetime + tlen_td / 3) or \
           np.max(rh_df['datenum'].values) < date2num(tdatetime + 2*tlen_td / 3) or len(rh_df.index) < 2:
            print('not enough data - continue')
            continue

        if smoothqc:
            presm = len(rh_df.index)
            #print(str(presm) + ' points pre smooth')
            idfw = np.where(rh_df['tane_dedt'].values > 0)[0]
            idbk = np.where(rh_df['tane_dedt'].values < 0)[0]
            smoothl = 5
            fwarr = rh_df.iloc[idfw]['refl_hgt'].astype('float')
            bkarr = rh_df.iloc[idbk]['refl_hgt'].astype('float')
            rhsmoothfw = mv_avg(np.array(fwarr), smoothl)
            rhsmoothbk = mv_avg(np.array(bkarr), smoothl)
            difffw = np.abs(rhsmoothfw - rh_df.iloc[idfw]['refl_hgt'].values)
            diffbk = np.abs(rhsmoothbk - rh_df.iloc[idbk]['refl_hgt'].values)
            filfw = difffw < 3 * np.std(difffw)
            filbk = diffbk < 3 * np.std(diffbk)
            #print('standard deviations are ' + str(np.std(difffw)) + ' and ' + str(np.std(diffbk)))
            idfw = idfw[filfw]
            idbk = idbk[filbk]
            ids_all = np.append(idfw, idbk)
            rh_df = rh_df.iloc[ids_all]
            #print('got rid of ' + str(presm - len(rh_df.index)) + ' points')

        if largetides:
            temp_dn = np.sort(rh_df['datenum'].values)
        else:
            temp_dn = np.sort(snrdt_df['datenum'].values)
        tfilter = np.logical_and(temp_dn > date2num(tdatetime + tlen_td / 3),
                                 temp_dn < date2num(tdatetime + 2 * tlen_td / 3))
        tfilter = np.where(tfilter)[0]
        try:
            if tfilter[0] != 0:
                tfilter = np.append(tfilter[0] - 1, tfilter)
            if tfilter[-1] != len(temp_dn) - 1:
                tfilter = np.append(tfilter, tfilter[-1] + 1)
            temp_dn = temp_dn[tfilter]
        except IndexError:
            print('not enough data - continue')
            continue
        maxtgap = np.max(np.ediff1d(temp_dn))
        mintgap = np.min(np.ediff1d(temp_dn))
        if mintgap < 0:
            print('issue - values not in order')
            continue
        print('max gap is ' + str(int(maxtgap * 24 * 60)) + ' minutes')

        if maxtgap > kspac:
            print('gap in data bigger than node spacing')
            # print('continue with risk of instabilities')
            continue

        knots = np.hstack((tdatenum * np.ones(bspline_order),
                          np.linspace(tdatenum, tdatenum_end, int(tlen / kspac + 1)),
                          tdatenum_end * np.ones(bspline_order)))

        def residuals_spectral_ls(sfacs):
            residuals = residuals_bsp_spectral(sfacs, bspline_order, knots, rh_df)
            return residuals

        sfacs_0 = np.nanmean(rh_df['refl_hgt'].values) * np.ones(int(tlen / kspac + bspline_order))
        # bounds = Bounds(rhlims[0], rhlims[1])
        ls_spectral = least_squares(residuals_spectral_ls, sfacs_0, method='trf', bounds=rhlims)
        # ls_spectral = least_squares(residuals_spectral_ls, sfacs_0, method='lm')
        sfacs_spectral = ls_spectral.x
        invout['sfacs_spectral'] = sfacs_spectral
        if snrfit:
            print('now doing Joakim Strandberg inversion')
            if largetides:
                sfacs_0 = sfacs_spectral
            else:
                sfacs_0 = np.median(rh_df['refl_hgt'])*np.ones(int(tlen / kspac + bspline_order))
            consts = len(satconsts)
            sfacs_0 = np.append(sfacs_0, np.zeros(consts * 2))
            if 'rough_in' in kwargs:
                sfacs_0 = np.append(sfacs_0, kwargs.get('rough_in'))

            def residuals_js_ls(sfacs):
                residuals = residuals_bsp_js(sfacs, bspline_order, knots, satconsts, signal, snrdt_df)
                return residuals

            ls_js = least_squares(residuals_js_ls, sfacs_0, method='lm')
            invout_js = ls_js.x
            sfacs_js = invout_js[:int(tlen / kspac + bspline_order)]
            params_js = invout_js[int(tlen / kspac + bspline_order):]
            invout['sfacs_js'] = sfacs_js
            invout['params_js'] = params_js
        tfilter = np.logical_and(rh_df['datenum'].values > date2num(tdatetime + tlen_td / 3),
                                 rh_df['datenum'].values < date2num(tdatetime + 2 * tlen_td / 3))
        rh_df = rh_df[tfilter]
        # print(str(len(rh_df.index)) + ' points')
        # saving kspac and tlen so that the invplot function doesn't require them as inputs
        invout['kspac'] = kspac
        invout['tlen'] = tlen
        invout['bspline_order'] = bspline_order
        Path(invdir).mkdir(parents=True, exist_ok=True)
        invfilestr = invdir + '/' + str(curdt.strftime("%y_%m_%d_%H_%M")) + '.pkl'
        f = open(invfilestr, 'wb')
        pickle.dump(invout, f)
        pickle.dump(rh_df, f)
        f.close()
        print('dumped a pickle')
    return


def invsnr_plot(sdatetime, edatetime, invdir, plotl=1/(24*10), cubspl=False, **kwargs):
    """
    plot output from the 'invsnr' function above and compare with tide gauge, if given as input
    the figure is saved to a file, 'reflh_test.png'
    :param sdatetime: start date and time, in datetime format
    :param edatetime: end date and time, in datetime format
    :param invdir: str to the directory(ies) with output from the invsnr function saved as pickle files
    IMPORTANT: you can take an average of the output from co-located antennas by using more than one string as an input
    if you just use one antenna then it must be contained in square brackets, e.g., ['path/to/inv/data/']
    :param plotl: the output temporal resolution of the bspline (in days, e.g., 60/86400 = 1 minute)
    :param kwargs:
    tgstring: string path to tide gauge data as a 'pickle' file, with variable tgdata,
    tgdata format is first column = datenum format, second column = sea level (m)
    see the function below 'read_canadian_tg'
    :return:
    """
    # first search for a file in directory on start day
    tfiles = listdir(invdir[0])
    tfiles_day = [tfile for tfile in tfiles if tfile[0:8] == str(sdatetime.strftime("%y_%m_%d"))]
    if len(tfiles_day) == 0:
        print('you need to choose the start date better')
        return
    startfilestr = invdir[0] + '/' + tfiles_day[0]
    # startfilestr = invdir[0] + '/' + str(sdatetime.strftime("%y_%m_%d_%H_%M")) + '.pkl'
    f = open(startfilestr, 'rb')
    invout = pickle.load(f)
    f.close()
    tlen = invout['tlen']
    kspac = invout['kspac']
    bspline_order = invout['bspline_order']
    tlen_td = datetime.timedelta(hours=int(tlen * 24))
    dispmissedmsg = True
    cntdir = -1
    for invdt in invdir:
        sfacs_spectral = np.empty(0)
        sfacs_js = np.empty(0)
        rh_df = pd.DataFrame({})
        tdatetime = sdatetime - tlen_td / 3
        while tdatetime < edatetime - tlen_td/3:
            tdatetime = tdatetime + tlen_td/3
            tfilestr = invdt + '/' + str(tdatetime.strftime("%y_%m_%d_%H_%M")) + '.pkl'
            inds = int(tlen / (3 * kspac) + 1)
            inde = int((2 * tlen) / (3 * kspac) + 1)
            if tdatetime == sdatetime:
                inds = 0
            elif tdatetime == edatetime - tlen_td/3:
                inde = int(tlen/kspac + bspline_order)
            if cubspl:
                inds = int(tlen / (3 * kspac))
                inde = int((2 * tlen) / (3 * kspac))
                if tdatetime == edatetime - tlen_td / 3:
                    inde = int((2 * tlen) / (3 * kspac)) + 1
                if tdatetime == sdatetime:
                    inds = int(tlen / (3 * kspac)) - 1
            try:
                f = open(tfilestr, 'rb')
                invout = pickle.load(f)
                if invout['tlen'] != tlen or invout['kspac'] != kspac:
                    print('the time window or knot spacing dont match - stopping')
                    return
                rh_dft = pickle.load(f)
                f.close()
                frames = [rh_df, rh_dft]
                rh_df = pd.concat(frames, ignore_index=True)
                # now get the scaling factors if they exist
                sfacs_spectralt = invout['sfacs_spectral'][inds:inde]
                if 'sfacs_js' in invout:
                    sfacs_jst = invout['sfacs_js'][inds:inde]
            except IOError:
                if dispmissedmsg:
                    print('missing data on ' + str(tdatetime.date()) + ' ' + str(tdatetime.time()) + ' putting nans')
                    dispmissedmsg = False
                sfacs_spectralt = np.empty(inde-inds)
                sfacs_spectralt[:] = np.nan
                if 'sfacs_js' in invout:
                    sfacs_jst = np.empty(inde-inds)
                    sfacs_jst[:] = np.nan
            sfacs_spectral = np.append(sfacs_spectral, sfacs_spectralt)
            if 'sfacs_js' in invout:
                sfacs_js = np.append(sfacs_js, sfacs_jst)
        if len(invdir) > 1:
            cntdir = cntdir + 1
            if cntdir == 0:
                sfacs_spectral_toavg = sfacs_spectral
                rh_df_toavg = rh_df
                if 'sfacs_js' in invout:
                    sfacs_js_toavg = sfacs_js
            else:
                sfacs_spectral_toavg = np.vstack((sfacs_spectral_toavg, sfacs_spectral))
                if 'sfacs_js' in invout:
                    sfacs_js_toavg = np.vstack((sfacs_js_toavg, sfacs_js))
                rh_df['antenna_id'] = cntdir
                frames = [rh_df_toavg, rh_df]
                rh_df_toavg = pd.concat(frames, ignore_index=True)
    if len(invdir) > 1:
        # DO THE AVERAGING HERE
        # but make it so that it would work the same with just one antenna
        meanhgts_spectral = np.empty((len(invdir)))
        if 'sfacs_js' in invout:
            meanhgts_js = meanhgts_spectral
        for cnt in range(0, cntdir+1):
            meanhgts_spectral[cnt] = np.nanmean(sfacs_spectral_toavg[cnt, :])
            if 'sfacs_js' in invout:
                meanhgts_js[cnt] = np.nanmean(sfacs_js_toavg[cnt, :])

            sfacs_spectral_toavg[cnt, :] = sfacs_spectral_toavg[cnt, :] + meanhgts_spectral[0] - \
                                           meanhgts_spectral[cnt]
            # now adjust the arc estimates, for plotting I guess
            rh_df_toavg.loc[rh_df_toavg['antenna_id'] == cnt, 'refl_hgt'] = \
                rh_df_toavg.loc[rh_df_toavg['antenna_id'] == cnt, 'refl_hgt'] + meanhgts_spectral[0] - \
                meanhgts_spectral[cnt]
            if 'sfacs_js' in invout:
                sfacs_js_toavg[cnt, :] = sfacs_js_toavg[cnt, :] + meanhgts_js[0] - meanhgts_js[cnt]
        rh_df = rh_df_toavg
        sfacs_spectral = np.nanmedian(sfacs_spectral_toavg, axis=0)
        if 'sfacs_js' in invout:
            sfacs_js = np.nanmedian(sfacs_js_toavg, axis=0)
    sdatenum = date2num(sdatetime)
    edatenum = date2num(edatetime)
    tplot = np.linspace(sdatenum, edatenum, int((edatenum - sdatenum) / plotl + 1))
    if not cubspl:
        knots = np.hstack((sdatenum - tlen / 3 * np.ones(bspline_order),
                           np.linspace(sdatenum - tlen / 3, edatenum + tlen / 3,
                                       int((edatenum - sdatenum + 2 * tlen / 3) / kspac + 1)),
                           edatenum + tlen / 3 * np.ones(bspline_order)))
        rh_spectral = interpolate.splev(tplot, (knots, sfacs_spectral, bspline_order))
    else:
        knots = np.linspace(sdatenum - kspac/2, edatenum + kspac/2, int((edatenum - sdatenum) / kspac + 2))
        tfilter = np.isnan(sfacs_spectral) == 0
        sfacs_spectral = sfacs_spectral[tfilter]
        knots = knots[tfilter]
        tfilter = np.logical_and(tplot[:] >= np.min(knots), tplot[:] <= np.max(knots))
        tplot = tplot[tfilter]
        cubspl_f = interpolate.interp1d(knots, sfacs_spectral, kind='cubic')
        rh_spectral = cubspl_f(tplot)

    fig, ax = plt.subplots(figsize=(8, 4))
    if 'tgstr' in kwargs:
        rh_df['refl_hgt'] = np.mean(rh_df['refl_hgt'].values) - rh_df['refl_hgt'].values
        tgstr = kwargs.get('tgstr')
        f = open(tgstr, 'rb')
        tgdata = pickle.load(f)
        tfilter = np.logical_and(tgdata[:, 0] >= sdatenum, tgdata[:, 0] < edatenum)
        tgdata = tgdata[tfilter, :]
        tgdata[:, 1] = tgdata[:, 1] - np.nanmean(tgdata[:, 1])
        ptg, = plt.plot_date(tgdata[:, 0], tgdata[:, 1], '-')
        ptg.set_label('tide gauge')
        # now calculate rms
        rh_spectral_rms = interpolate.splev(tgdata[:, 0], (knots, sfacs_spectral, bspline_order))
        if cubspl:
            tfilter = np.logical_and(tgdata[:, 0] >= np.min(knots), tgdata[:, 0] <= np.max(knots))
            rh_spectral_rms = cubspl_f(tgdata[tfilter, 0])
        tf = ~np.isnan(rh_spectral_rms)
        rh_spectral_rms = np.nanmean(rh_spectral_rms) - rh_spectral_rms
        rh_spectral = np.nanmean(rh_spectral) - rh_spectral
        if not cubspl:
            rms_spectral = np.sqrt(np.mean((rh_spectral_rms[tf] - tgdata[tf, 1]) ** 2))*100
        else:
            rms_spectral = np.sqrt(np.mean((rh_spectral_rms[tf] - tgdata[tfilter, 1]) ** 2)) * 100
        print('rms_spectral is %.3f cm' % rms_spectral)
        # corr_spectral = pearsonr(rh_spectral_rms[tf], tgdata[tf, 1])
        # print('corr_spectral is %.3f' % corr_spectral[0])
    parc, = plt.plot_date(rh_df['datenum'], rh_df['refl_hgt'], '.')
    parc.set_label('arcs')
    pspec, = plt.plot_date(tplot, rh_spectral, '-')
    pspec.set_label('b-spline fit of arcs')
    if 'sfacs_js' in invout:
        rh_js = interpolate.splev(tplot, (knots, sfacs_js, bspline_order))
        if 'tgstr' in kwargs:
            rh_js_rms = interpolate.splev(tgdata[:, 0], (knots, sfacs_js, bspline_order))
            rh_js_rms = np.nanmean(rh_js_rms) - rh_js_rms
            rh_js = np.nanmean(rh_js) - rh_js
            tf = ~np.isnan(rh_js_rms)
            rms_js = np.sqrt(np.mean((rh_js_rms[tf] - tgdata[tf, 1]) ** 2))*100
            print('rms_js is %.3f cm' % rms_js)
            corr_js = pearsonr(rh_js_rms[tf], tgdata[tf, 1])
            # print('corr_js is %.3f' % corr_js[0])
        pjs, = plt.plot_date(tplot, rh_js, '-')
        pjs.set_label('inverse modelling of SNR')
    print(str(round(len(rh_df.index) / (edatenum - sdatenum))) + ' points per day')
    dformat = DateFormatter('%Y-%m-%d %H:%M')
    ax.xaxis.set_major_formatter(dformat)
    ax.set_xlim(sdatenum, edatenum)
    ax.set_xticks(np.linspace(sdatenum, edatenum, int((edatenum-sdatenum)*2+1)))
    ax.legend()
    plt.savefig('reflh_test.png')
    plt.close()
    return rh_df


def read_canadian_tg(tgstr, tgstr_out):
    """
    read canadian tide gauge data, saved from
    http://www.isdm-gdsi.gc.ca/isdm-gdsi/twl-mne/maps-cartes/inventory-inventaire-eng.asp?user=isdm-gdsi&region=MEDS&tst=1&perm=0
    :param tgstr: string to tide gauge data, csv file
    :param tgstr_out: string for pickle file to save the data
    :return:
    """
    # note that I had to add a comma after all the header rows to make it work
    tgcsv = pd.read_csv(tgstr, header=7, names=['date', 'slvl', 'nan'])
    tgslvl = np.array(tgcsv['slvl'].values, dtype=float)
    dates = tgcsv['date']
    datet = [datetime.datetime.strptime(date, '%Y/%m/%d %H:%M') for date in dates]
    tgdatenum = date2num(datet)
    tgdata = np.column_stack((tgdatenum, tgslvl))
    f = open(tgstr_out, 'wb')
    pickle.dump(tgdata, f)
    f.close()
    return


def mv_avg(arr, smoothl):
    if (smoothl % 2) == 0:
        print('smoothl must be odd number')
        return
    smoothout = np.empty((len(arr)))
    cumarr = np.cumsum(arr)
    for i in range(0, len(arr)):
        sti = i - int((smoothl - 1) / 2) - 1
        eni = i + int((smoothl - 1) / 2)
        if sti < 0:
            cumsums = 0
            sti = -1
        else:
            cumsums = cumarr[sti]
        if eni > len(arr) - 1:
            eni = len(arr) - 1
        smoothout[i] = (cumarr[eni] - cumsums) / (eni - sti)
    return smoothout




