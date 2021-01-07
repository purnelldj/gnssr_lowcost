addpath('functions')

clear all

outdir='GRE_1050_10h_06h_05s_obs_19oct_norm_80220_lt';
stations=cellstr({'rv3sa','rv3sb','rv3sc','rv3sd'});
startdate=datenum(2020,9,11);
enddate=datenum(2020,10,9);
makefig=1;
roughnessplot=0;
belly=0;
if belly==0
    datastr='~/data/';
else
    datastr='/home/dpurnell/projects/ctb-ng50/dpurnell/data/';
end
%tgstring='';
%tgstring=[datastr,'sab2/tg_jan1_jun15.mat'];
tgstring=[datastr,'rv3s/tg_sep9_oct10.mat'];
%tgstring=[datastr,'rv3s/tg_aug6_aug13.mat'];
%tgstring=[datastr,'nypm/nypmtg_sept2020.mat'];

if strcmp(outdir(7),'d')
kspac=str2double(outdir(9:10))/10/24;
else
kspac=str2double(outdir(10:11))/10/24;
end
if strcmp(outdir(16),'h')
tlen=str2double(outdir(14:15))/24;
else
tlen=3;
end
plotl=15/(24*60);
for jj=1:numel(stations)
invdir(jj)=cellstr({[datastr,char(stations(jj)),'/',outdir]});
end

[t_rh,rh_invjs,rh_invpre,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,makefig,roughnessplot);



