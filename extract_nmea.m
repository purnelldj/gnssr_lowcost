%clear all;

for ii=1:4
station='rv3s';
antennaid=upper(char(96+ii));
filenums=9:11;
pwdstr=pwd;
snrdir=['data/',station,lower(antennaid),'/snr/'];
sp3dir='~/data/sp3';
elv_lims=[0 90];
azi_lims=[70 240];
splice=15;
sp3option=3;

addpath('functions')
for ff=1:numel(filenums)
    disp(['doing file number ',num2str(filenums(ff))])
    nmeafile=['~/data/',station,'/rawdata_10oct/gpslog',antennaid,...
        num2str(filenums(ff))];
    nmea2snr_fun(nmeafile,sp3dir,sp3option,elv_lims,azi_lims,splice,snrdir)
end

end


