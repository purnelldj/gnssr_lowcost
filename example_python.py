import lc_fun
import datetime

# all code written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost
# run using python v3.6

# 1. first extract SNR data from NMEA files for analysis
# more precise elevation and azimuth angle data is extracted from orbit '.sp3' files

nmeastr = 'data/rv3s/rawdata_github/gpslogA_sample'
sp3dir = 'data/sp3/COD'
snrdir = 'data/rv3s/a/snr'
lc_fun.nmea2snr(nmeastr, snrdir, sp3dir=sp3dir, tempres=15)

# 2. now analyze the SNR data

snrdir = 'data/rv3s/a/snr'
invdir = 'data/rv3s/a/invout'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
kspac = 1/24
tlen = 6/24
rhlims = [3.5, 6]
lc_fun.invsnr(sdatetime, edatetime, snrdir, invdir, kspac, tlen, rhlims, pktnlim=15, smoothqc=True,
              elvlims=[10, 50], azilims=[80, 220], tempres=15, rough_in=0.001, snrfit=True, arctlim=1200,
              normalize=True)

# 3. then plot the output from inverse analysis

# put square brackets around even if just one entry
invdir = ['data/rv3s/a/invout']  # SQUARE BRACKETS ARE NECESSARY:
# can plot more than one antenna by separating strings with comma
tgstr = 'data/rv3s/rv3sep.pkl'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
lc_fun.invsnr_plot(sdatetime, edatetime, invdir, tgstr=tgstr)

# should be a figure output rel_hgt.png

