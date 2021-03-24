import lc_fun
import datetime

# 1. first extract SNR data from NMEA files for analysis
# more precise elevation and azimuth angle data is extracted from orbit '.sp3' files

nmeastr = 'data/rv3s/rawdata_github/gpslogD_sample'
sp3dir = 'data/sp3/COD'
snrdir = 'data/rv3s/d/snr'

lc_fun.nmea2snr(nmeastr, sp3dir, snrdir, tempres=5/86400)

exit()

# 2. now analyze the SNR data

snrdir = 'data/rv3s/d/snr'
invdir = 'data/rv3s/d/invout'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
kspac = 1/24
tlen = 6/24
rhlims = [4, 6]
lc_fun.invsnr(sdatetime, edatetime, snrdir, invdir, kspac, tlen, rhlims, pktnlim=15, arctlim=1200/86400,
              elvlims=[10, 40], azilims=[80, 220], tempres=15, rough_in=0.001, snrfit=True, normalize=True)

exit()

# 3. then plot the output from inverse analysis

# put square brackets around even if just one entry
invdir = ['data/rv3s/a/invout', 'data/rv3s/b/invout', 'data/rv3s/c/invout', 'data/rv3s/d/invout']
tgstr = 'data/rv3s/rv3sep.pkl'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
lc_fun.invsnr_plot(sdatetime, edatetime, invdir, tgstr=tgstr)

exit()
