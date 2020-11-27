% to analyze a period of SNR data for a certain site
clear all;

station='rv3sa';
startdate=datenum(2020,5,17);
enddate=datenum(2020,5,25); % including the last day
satconsts=[1,1,1];
outdir='obs_stats_3060';
snrdir='snr_90';
tropd=0;
splice=15;
tempsnr=0;
templsp=0;
tempfresnel=0;
datastr='~/data/';
tdatenum=startdate-1;
addpath('functions')

% LOAD SNR DATA
while tdatenum<enddate
tdatenum=tdatenum+1;
curdt=datetime(tdatenum,'convertfrom','datenum');
disp(char(curdt))
strday=char(datetime(curdt,'format','DDD'));
stryr=char(datetime(curdt,'format','yy'));
if exist([datastr,station,'/',snrdir,'/',num2str(tdatenum),'.mat'])~=0
    load([datastr,station,'/',snrdir,'/',num2str(tdatenum),'.mat'])
elseif exist([datastr,station,'/snr/',station,strday,'0.',stryr,'snr.mat'])~=0
    load([datastr,station,'/snr/',station,strday,'0.',stryr,'snr.mat'])
    snr_data=snrdata;
elseif exist([datastr,station,'/snr/',station,strday,'0.',stryr,'snr'])~=0
    snr_data=dlmread([datastr,station,'/snr/',station,strday,'0.',stryr,'snr']);
else
    disp('snr data does not exist!')
    continue
end
%%%%%%%%%%%%%%%
[slvlr,lspy] = analyzesnr_fun(station,snr_data,tdatenum,satconsts,...
    tropd,splice,tempsnr,templsp,tempfresnel);
% now save the output
if exist([datastr,station,'/',outdir])==0
    mkdir([datastr,station,'/',outdir])
end
save([datastr,station,'/',outdir,'/',num2str(tdatenum),'.mat'],'slvlr','lspy')
clear slvlr lspy
end

