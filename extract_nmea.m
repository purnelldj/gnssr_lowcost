%clear all;

for ii=1:4
station='rv3s';
antennaid=upper(char(96+ii));
filenums=1:11;
pwdstr=pwd;
snrdir=['data/',station,lower(antennaid),'/snr/'];
sp3dir='~/data/sp3';
elvlims=[0 90];
azilims=[70 240];
decimate=15;
sp3option=3;

addpath('functions')
for ff=1:numel(filenums)
    disp(['doing file number ',num2str(filenums(ff))])
    nmeafile=['~/data/',station,'/rawdata_12nov/gpslog',antennaid,...
        num2str(filenums(ff))];
    nmea2snr(nmeafile,sp3dir,sp3option,elvlims,azilims,decimate,snrdir);
end

end


