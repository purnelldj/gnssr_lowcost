% to obtain reflector height time series via inverse modelling
% as per Strandberg et al. (2016) (invsnr.m)

clear all;

%station='rv3sa';
station='rotg';
disp(':)(:')
disp(station)
disp(':)(:')
outdir='invout_test';
startdate=datenum(2015,1,1);
enddate=datenum(2015,1,5); % includes last day

datastr='~/data/';
snrdir=cellstr({['~/data/',station,'/snr']});

% basic station parameters
% rv3sa
%staxyz=[1323539.9905,-4207749.0481,4591443.2485];
%elvlims=[10 40];
%azilims=[80 220];
%rhlims=[2 6];

% rotg
staxyz=[4205957.7325,-291580.9908,4770001.1825];
azilims=[150 270];
elvlims=[20 30];
rhlims=[2 14];

satconsts=[1 1 0]; % GPS, GLO, GAL; 1 or 0
arclims=[300 inf]; % upper / lower time limit of arcs
kspac=1.5/24; % in days
tlen=9/24; % in days
dt=10; % in seconds
sig=1; % l1 or l2
roughin=0.01; % in meters
largetides=1;
skipjs=1;

addpath('functions')

tdatenum=startdate-tlen/3;
while round(tdatenum,10,'significant')<round(enddate+1-tlen/3,10,'significant')
tdatenum=tdatenum+tlen/3;
[sfacsjs,sfacspre,xinit,hinit,consts_out,roughout] = invsnr(tdatenum,snrdir,staxyz,...
    elvlims,azilims,rhlims,kspac,tlen,dt,satconsts,sig,arclims,largetides,roughin,skipjs);
if exist([datastr,station,'/',outdir])==0
    mkdir([datastr,station,'/',outdir]);
end
if numel(hinit)==1 && isnan(hinit)
    disp('not saving')
else
save([datastr,station,'/',outdir,'/',num2str(round(tdatenum,10,'significant')),'.mat'],...
    'sfacsjs','sfacspre','xinit','hinit','consts_out','roughout')
end

end

return

%% now to plot bspline

addpath('functions')

clear all
close all

datastr='~/data/';
invdir=cellstr({[datastr,'/rotg/GR0_2030_15h_09h_10s_obs_nojs_150270_l2'],[datastr,'/rotg/GR0_2030_15h_09h_10s_obs_nojs_150270']});
startdate=datenum(2015,1,1);
enddate=datenum(2015,1,5);
tgstring='~/data/rotg/tg_2015_10min.mat';

kspac=1.5/24;
tlen=9/24;
plotl=10/(24*60);
rhlims=[2 14];
plotrh=1;
plotvdc=0;

[t_rh,rh_invjs,rh_invpre,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,rhlims,plotrh,plotvdc);

