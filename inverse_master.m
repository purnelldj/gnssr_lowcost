% to obtain reflector height time series via inverse modelling
% as per Strandberg et al. (2016) (invsnr.m)

clear all;

%%%%%%%%%%%%%%%%%%%%%%%
% CHANGE ONLY HERE
%for ii=2:4
belly=0;
station=({'rv3sa','rv3sb','rv3sc','rv3sd'});
disp(':)(:')
disp(station)
disp(':)(:')
outdir='GRE_2050_10h_72h_15s_obs_19oct_norm_multi_abcd';
startdate=datenum(2020,9,18);
enddate=datenum(2020,9,21);
%%%%%%%%%%%%%%%%%%%%%%%

if belly==1
    datastr='/home/dpurnell/projects/ctb-ng50/dpurnell/data/';
    addpath('/home/dpurnell/projects/ctb-ng50/dpurnell/work/gnssr_matlab/functions')
    cd '/home/dpurnell/projects/ctb-ng50/dpurnell/work/gnssr_matlab'
else
    datastr='~/data/';
    addpath('/Users/Dave/Library/Mobile Documents/com~apple~CloudDocs/work/gnssr_matlab/functions')
    cd '/Users/Dave/Library/Mobile Documents/com~apple~CloudDocs/work/gnssr_matlab'
end
% these automatically change with the above
if strcmp(outdir(22:24),'obs')
    if size(station,2)<1
        snrdir=cellstr({[datastr,station,'/snr']});
    else
        for ii=1:size(station,2)
            snrdir(ii)=cellstr({[datastr,char(station(ii)),'/snr']});
        end
    end
else
snrdir=cellstr({[datastr,station,'/synth_gpsglo']});
end
satconsts=[0 0 0];
if strcmp(outdir(1),'G')==1
    satconsts(1)=1;
end
if strcmp(outdir(2),'R')==1
    satconsts(2)=1;
end
if strcmp(outdir(3),'E')==1
    satconsts(3)=1;
end
altelvlims=[str2double(outdir(5:6)) str2double(outdir(7:8))];
kspac=str2double(outdir(10:11))/10/24;
decimate=str2double(outdir(18:19)); % in seconds
tlen=str2double(outdir(14:15))/24;
roughin=0.001; % 0.001 for low-cost
largetides=0;
skipjs=1;

tdatenum=startdate-tlen/3;
while round(tdatenum,10,'significant')<round(enddate+1-tlen/3,10,'significant')
tdatenum=tdatenum+tlen/3;
%return
[sfacsjs,sfacspre,hinit,xinit,consts_out,roughout] = invsnr_norm_multi(tdatenum,station,snrdir,...
    kspac,tlen,decimate,satconsts,altelvlims,largetides,roughin,skipjs);
if exist([datastr,char(station(1)),'/',outdir])==0
    mkdir([datastr,char(station(1)),'/',outdir]);
end
if numel(hinit)==1 && isnan(hinit)
    disp('not saving')
else
save([datastr,char(station(1)),'/',outdir,'/',num2str(round(tdatenum,10,'significant')),'.mat'],...
    'sfacsjs','sfacspre','hinit','xinit','consts_out','roughout')

end

clearvars -except jj

%end
