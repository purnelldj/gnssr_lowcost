% to perform basic analysis (spectral analysis) on snr arcs

clear all
addpath('matlab_functions')

station='rv3s/a';
startdate=datenum(2020,9,13);
enddate=datenum(2020,9,13);
snrdir=['data/',station,'/snr'];
outdir=['data/',station,'/arcs_stats'];

% then choose these parameters
elvlims = [10 50]; % lower / upper limit
azilims = [80 220]; % lower / upper limit
rhlims = [3.5 6]; % ignore reflector heights smaller or larger than these limits
satconsts=[1 1 1]; % GPS, GLO, GAL; 1 or 0
staxyz = [1323539.0504, -4207748.7536, 4591443.7857];
dt=15;
sig=1;
normalize=0;

pktnlim=15;
arclims=1200;

snrfigs=1;
lspfigs=1;
gfresnel=0;

tdatenum=startdate-1;
while tdatenum<enddate
tdatenum=tdatenum+1;
curdt=datetime(tdatenum','convertfrom','datenum');
disp(char(curdt))
load([snrdir,'/',num2str(tdatenum),'.mat'])
if size(snr_data,2)==8
    snr_data(:,9)=tdatenum+snr_data(:,4)./86400;
end

[rh_stats,~] = snr2arcs(snr_data,staxyz,elvlims,azilims,...
    rhlims,dt,satconsts,sig,arclims,normalize,pktnlim,snrfigs,lspfigs,gfresnel);
if ~exist(outdir,'dir')
    mkdir(outdir);
end
if numel(rh_stats)==1 && isnan(rh_stats)
    disp('not saving')
else
save([outdir,'/',num2str(round(tdatenum,10,'significant')),'.mat'],'rh_stats')
end

end
rh_stats = sortrows(rh_stats, 3);


return
figure('visible', 'on')
scatter(rh_stats(:,1), rh_stats(:,2))
