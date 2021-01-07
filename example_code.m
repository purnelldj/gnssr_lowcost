% here are some examples for how to run various codes
% PLEASE NOTE the descriptions of inputs are given in the function descriptions
% ALL WRITTEN BY DAVID PURNELL copyright 2020
% to go with the article:
% "Quantify the Uncertainty in Ground-Based GNSS-Reflectometry Sea Level
% Measurements" (2020) (IEEE JSTARS) by Purnell et al.

% PLEASE SEE THE README FILE FIRST AND FOLLOW INSTRUCTIONS
% THESE CODES WONT WORK UNLESS YOU LOOK AT THE README FILE

addpath('functions')

%% first extract snr data from rinex file (rinex2snrfile.m)

clear all;

station='sc02';
startdate=datenum(2015,1,1);
enddate=datenum(2015,1,5); % including the last day
outdir='snr'; % output directory
sp3option=1;
elv_lims=[0 40];
azi_lims=[0 360];
dt=15; % the frequency of data (in seconds)

% note this is just to get staxyz, lat and lon
run(['functions/station_inputs/',station,'_input.m'])

addpath('functions')
tdatenum=startdate-1;
datastr='data/';
% LOAD SNR DATA
while tdatenum<enddate
tdatenum=tdatenum+1;
curdt=datetime(tdatenum,'convertfrom','datenum');
curjd=juliandate(curdt);
[gpsw,sow,~]=jd2gps(curjd);
dow=sow/86400;
disp(char(curdt))
strday=char(datetime(curdt,'format','DDD'));
stryr=char(datetime(curdt,'format','yy'));
stryrl=char(datetime(curdt,'format','yyyy'));

if exist([datastr,station,'/obs/',station,strday,'0.',stryr,'o'],'file')==2
    obsstr=[datastr,station,'/obs/',station,strday,'0.',stryr,'o'];
else
    disp('snr data does not exist!')
    continue
end

% choose an option for the sp3 files
% option 3 best for multi-GNSS
% get sp3 files using functions/getfiles/getsp3.m
if sp3option==1
% OPTION 1
sp3str=[datastr,'sp3/com',num2str(gpsw),num2str(round(dow)),'.sp3'];
elseif sp3option==2
% OPTION 2
sp3str=[datastr,'sp3/igs',num2str(gpsw),num2str(round(dow)),'.sp3'];
elseif sp3option==3
% OPTION 3
sp3str=[datastr,'sp3/COD/COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
elseif sp3option==4
% OPTION 4
sp3str=[datastr,'sp3/GFZ/GFZ0MGXRAP_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
end

[snr_data]=rinex2snrfile(obsstr,sp3str,elv_lims,azi_lims,staxyz,plat,plon);

%%%%%%%%%%%%%%%
% now save the output
if exist([datastr,station,'/',outdir])==0
    mkdir([datastr,station,'/',outdir])
end

% saving SNR files here
save([datastr,station,'/',outdir,'/',num2str(tdatenum),'.mat'],'snr_data')
clear slvlr lspy
end

%% then analyze some observed SNR data (analyzesnr_fun.m)

station='sc02';
satconsts=[1,0,0]; % just GPS
tropd=0;
decimate=15;
tempsnr=0;
templsp=0;
tempfresnel=0;
startdate=datenum(2015,1,1);
enddate=datenum(2015,1,5); % including last day
outdir='obs_stats'; % output directory

slvlr_all=[];
tdatenum=startdate-1;
datastr='data/';
while tdatenum<enddate
tdatenum=tdatenum+1;
% LOAD SNR DATA
curdt=datetime(tdatenum,'convertfrom','datenum');
disp(char(curdt))
strday=char(datetime(curdt,'format','DDD'));
stryr=char(datetime(curdt,'format','yy'));
load([datastr,station,'/snr/',num2str(tdatenum),'.mat'])
snrfile=snr_data;
%%%%%%%%%%%%%%%
[slvlr,lspy] = analyzesnr_fun(station,snrfile,tdatenum,satconsts,tropd,decimate,tempsnr,templsp,tempfresnel);
% now save the output

if exist([datastr,station,'/',outdir])==0
    mkdir([datastr,station,'/',outdir])
end

save([datastr,station,'/',outdir,'/',num2str(tdatenum),'.mat'],'slvlr','lspy')
slvlr_all=[slvlr_all;slvlr];
end

% plotting unadjusted reflector height estimates from spectral analysis
scatter(slvlr_all(:,1),slvlr_all(:,7),'b+')

%% to obtain reflector height time series via inverse modelling as per Strandberg et al. (2016) (invsnr.m)

clear all

startdate=datenum(2015,1,2);
enddate=datenum(2015,1,2);
station='sc02';
%%%%%%%%%%%%%
% observed
%snrdir=cellstr({'data/sc02/snr'});
%outdir='data/sc02/inv_test';
% synthetic
snrdir=cellstr({'data/sc02/synth_data'});
outdir='data/sc02/inv_test_synth';
%%%%%%%%%%%%%
kspac=3/24; % in days
tlen=18/24; % in days
decimate=0; % in seconds
satconsts=[1 0 0];
roughin=0.1;
altelvlims=[];
largetides=1;
skipjs=0;

tdatenum=startdate-tlen/3;
while round(tdatenum,10,'significant')<round(enddate+1-tlen/3,10,'significant')
tdatenum=tdatenum+tlen/3;
[sfacsjs,sfacspre,hinit,xinit,consts_out,roughout] = invsnr(tdatenum,station,snrdir,...
    kspac,tlen,decimate,satconsts,altelvlims,largetides,roughin,skipjs);
if exist(outdir)==0
    mkdir(outdir);
end
save([outdir,'/',num2str(round(tdatenum,10,'significant')),'.mat'],...
    'sfacsjs','sfacspre','hinit','xinit','consts_out','roughout')
end

%% to plot inverse reflector height estimates and compare with tide gauge or reference (invsnr_plot.m)

clear all

startdate=datenum(2015,1,2);
enddate=datenum(2015,1,2);
invdir=cellstr({'data/sc02/inv_test_synth'});
kspac=3/24;
tlen=18/24;
plotl=1/(24*60);
tgstring='data/sc02/tg_2015_6min.mat';
makefig=1;
roughnessplot=0;

[t_rh,rh_invjs,rh_invpre,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,makefig,roughnessplot);


