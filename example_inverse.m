% to obtain reflector height time series via inverse modelling
% as per Strandberg et al. (2016) (invsnr.m)

clear all;

station='rv3sa';
disp(':)(:')
disp(station)
disp(':)(:')
outdir='invout_test';
startdate=datenum(2020,10,9);
enddate=datenum(2020,10,10);

datastr='data/';
addpath('functions')

% basic station parameters
elvlims=[10 40];
azilims=[80 220];
staxyz=[1323539.9905,-4207749.0481,4591443.2485];
rhlims=[2 6];

snrdir=cellstr({['data/',station,'/snr']});
satconsts=[1 1 1];
arclims=[300 inf];
kspac=1/24; % in days
tlen=6/24; % in days
dt=15; % in seconds
roughin=0.001; % in meters
largetides=1;
skipjs=1;
sig=1; % l1 or l2

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
