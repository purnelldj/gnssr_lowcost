import lc_fun
import datetime
import pickle
import matplotlib.pyplot as plt

# 4. for looking at arcs using snr2arcs

snrdir = 'data/rv3s/a/snr'
curdt = datetime.datetime(2020, 9, 13)
rhlims = [3.5, 6]
snrfilestr = snrdir + '/' + str(curdt.strftime("%y_%m_%d")) + '.pkl'
f = open(snrfilestr, 'rb')
snr_df = pickle.load(f)
f.close()
rh_df, _ = lc_fun.snr2arcs(snr_df, rhlims, pktnlim=11.25, elvlims=[10, 50], azilims=[80, 220], arctlim=1200, tempres=15)
#fig, ax = plt.subplots(figsize=(8, 4))
#ax.plot(rh_df['datenum'], rh_df['refl_hgt'], '.')
print(len(rh_df.index))
#print(rh_df[rh_df['sat_prn'] == 'G10']['datenum'].values)

exit()

# 1. first extract SNR data from NMEA files for analysis
# more precise elevation and azimuth angle data is extracted from orbit '.sp3' files

nmeastr = 'data/rv3s/rawdata_github/gpslogD_sample'
sp3dir = 'data/sp3/COD'
snrdir = 'data/rv3s/d/snr'

lc_fun.nmea2snr(nmeastr, sp3dir, snrdir, tempres=15)

exit()

# 3. then plot the output from inverse analysis

# put square brackets around even if just one entry
invdir = ['data/rv3s/a/invout_nojs_tmp']
tgstr = 'data/rv3s/rv3sep.pkl'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
lc_fun.invsnr_plot(sdatetime, edatetime, invdir, tgstr=tgstr)

exit()

# 2. now analyze the SNR data

snrdir = 'data/rv3s/a/snr'
invdir = 'data/rv3s/a/invout_nojs_tmp'
sdatetime = datetime.datetime(2020, 9, 13)
edatetime = datetime.datetime(2020, 9, 14)
kspac = 1/24
tlen = 6/24
rhlims = [3.5, 6]
lc_fun.invsnr(sdatetime, edatetime, snrdir, invdir, kspac, tlen, rhlims, pktnlim=0, smoothqc=False,
              elvlims=[10, 50], azilims=[80, 220], tempres=15, rough_in=0.001, snrfit=False,
              normalize=True)

exit()
