% to analyze low-cost GNSS data and analyze for reflectometry
% using some example data from 'short' antenna array at Trois-Rivières

% first add directory that contains all the codes used below to path
addpath('functions')

%% analyze the raw output from low-cost antennas to extract SNR data
% example data is too large for github - skip to next part
% however will upload along with paper
% also see nmea2snr.m description

clear all

for ii = 1:1 % loop for antenna A, B, C and D
antennaid = upper(char(96+ii));
filenums = 1:1;
station = 'rv3s'; % name for station
snrdir = ['data/',station,'/',lower(antennaid),'/snr/']; % output directory
sp3dir = 'data/sp3'; % directory with sp3 orbit files
elvlims = [0 80]; % elevation limits
azilims = [80 220]; % azimuth limits
tempres = 15; % choose temporal resolution in seconds (raw data is 1 second)
sp3option = 3; % see in nmea2snr.m for different options

for ff = 1:numel(filenums)
    disp(['doing file number ',num2str(filenums(ff))])
    nmeafile = ['data/rv3s/gpslog',antennaid,num2str(filenums(ff))];
    nmea2snr(nmeafile,sp3dir,sp3option,elvlims,azilims,tempres,snrdir);
end

end

return

%% now extract some station parameters

clear all

statsdir='data/rv3s/a/snr';

files=dir(statsdir);

heights=[];
lats=[];
lons=[];
staxyzs=[];
names={files.name};
for ii=1:size(names,2)
    tmpfile=char(names(ii));
    if numel(tmpfile)>4
    if strcmp(tmpfile(end-3:end),'.mat')
        load([statsdir,'/',tmpfile])
        lats=[lats lat];
        lons=[lons lon];
        heights=[heights height];
        staxyzs=[staxyzs; staxyz];
    end
    end
end

format long
disp(['lat is ',num2str(nanmedian(lats))]);
disp(['lon is ',num2str(nanmedian(lons))]);
disp(['height is ',num2str(nanmedian(heights))]);
disp(['staxyz is ',num2str(nanmedian(staxyz,1))]);

%% do inverse modelling for one antenna at a time
% method adapted from Strandeberg et al. (2016)
% for more info see 'invsnr.m'

clear all

outdir = 'data/rv3s/a/invout_test';
snrdir = cellstr({'data/rv3s/a/snr'});
startdate = datenum(2020,9,11);
enddate = datenum(2020,9,12); % includes last day

% taken from above
% for rv3sa
staxyz = [1323539.5944,-4207748.858,4591444.0879];
% then choose these parameters
elvlims = [10 40]; % lower / upper limit
azilims = [80 220]; % lower / upper limit
rhlims = [4 6]; % ignore reflector heights smaller or larger than these limits
satconsts=[1 1 1]; % GPS, GLO, GAL; 1 or 0
arclims=[1200 inf]; % lower / upper time limit of arcs
kspac=1/24; % bspline knot spacing, in days
tlen=6/24; % bspline time window length, in days
dt=15; % in seconds
sig=1; % l1 or l2
roughin=0.001; % in meters
largetides=1; % if there are significant tides at site
skipjs=0; % 0 to do least squares fit as per Strandeberg et al. (2016)

addpath('functions')

tdatenum=startdate-tlen/3;
while round(tdatenum,10,'significant')<round(enddate+1-tlen/3,10,'significant')
tdatenum=tdatenum+tlen/3;
[sfacsjs,sfacspre,xinit,hinit,consts_out,roughout] = invsnr(tdatenum,snrdir,staxyz,...
    elvlims,azilims,rhlims,kspac,tlen,dt,satconsts,sig,arclims,largetides,roughin,skipjs);
if exist(outdir)==0
    mkdir(outdir);
end
if numel(hinit)==1 && isnan(hinit)
    disp('not saving')
else
save([outdir,'/',num2str(round(tdatenum,10,'significant')),'.mat'],...
    'sfacsjs','sfacspre','xinit','hinit','consts_out','roughout')
end

end

%% now to plot the output from inverse modelling

clear all; close all

invdir=cellstr({'data/rv3s/a/invout_test'}); % combine solutions by putting more than one directory
startdate=datenum(2020,9,11);
enddate=datenum(2020,9,12);
tgstring='data/rv3s/tg_sep9_oct10.mat';

kspac=1/24;
tlen=6/24;
plotl=10/(24*60);
rhlims=[4 6];
plotrh=1;
plotvdc=0;

[t_rh,rh_invjs,rh_invpre,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,rhlims,plotrh,plotvdc);

