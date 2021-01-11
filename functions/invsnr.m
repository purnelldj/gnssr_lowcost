function [sfacsjs,sfacspre,xinit,hinit,consts_out,roughout] = invsnr(tdatenum,snrdir,staxyz,...
    elvlims,azilims,rhlims,kspac,tlen,dt,satconsts,sig,arclims,largetides,roughin,skipjs)


%%
% this code is for inverse modelling of SNR data to get GNSS-R measurements
% as per Strandberg et al. (2016)

% TEST

% INPUTS
% tdatenum: day or time in datenum format
% snrdir: path to SNR data, should also be in cellstr format
% staxyz: approx xyz of statiom
% elvlims: elevation angle limits
% azilims: azimuth angle limits
% rhlims: reflector height limits min/max (m)
% kspac: average node spacing in days (e.g., 2/24 is 2 hours)
% tlen: length of window for analysis in days, split into 3 and the middle
% period is saved (e.g., setting to 3 means that the middle day is saved)
% dt: temporal resolution of input SNR data (in seconds), if larger than
% the actual resolution then will remove data to convert
% satconsts: 1 by 3 double, 1 or 0 to include satellite constellations
% GPS, GLONASSS, GALILEO (e.g., [1 0 1] for GPS and GALILEO)
% sig: 1 for l1 or 2 for l2
% arclims: 1 by 2 double, min and max time limits (in seconds) for a satellite arc to be
% included, suggested default [300 inf]
% altelvlims: to overwrite the elevation limits in the station input file
% and use alternate ones, e.g., [5 15] for 5 to 15 degrees
% largetides: changes the initial guess for node values (1 or 0)
% roughin: initial sfc roughness guess ('s') in metres for least squares
% adjustment, set to '' if don't want to use roughness
% skipjs: skip the least squares adjustment as per strandberg et. al and
% instead just fit the b-spline curve using spectral analysis estimates

% OUTPUTS
% these outputs are to go with the function 'invsnr_plot.m'
% sfacsjs: node values estimated from Strandberg et al. analysis
% sfacspre: node values estimated by using spectral analysis adjustment
% hinit: initial spectral analysis estimates
% xinit: timing of spectral analysis estimates
% consts_out: the other variables estimated as part of least squares
% adjustment
% roughout: the roughness parameter that is estimated in least squares
% adjustment

p=2; % bspline order
stdfac=3;

snrfigs=0;
lspfigs=0;
gfresnel=0;

gps=satconsts(1);
glo=satconsts(2);
gal=satconsts(3);

sfacsjs=NaN;
sfacspre=NaN;
xinit=NaN;
hinit=NaN;
consts_out=NaN;
roughout=NaN;

pwdstr=pwd;
addpath([pwdstr,'/functions/'])
addpath([pwdstr,'/functions/bspline'])

tdatenum=tdatenum-tlen/3;
 
curdt=datetime(tdatenum+tlen/3,'convertfrom','datenum');
disp(char(curdt))

% work out if need two days or just one
mlen=1;
if tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))>0
    mlen=tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))+1;
    if mod(tdatenum+tlen,1)==0
        mlen=mlen-1;
    end
end

% AT THIS POINT THE tdatenum AND CURJD ARE THE SAME

% try getting day in order of station first, then after detrended, organise
% by time, then that's it
snr_data_all=[];
for ll=1:numel(snrdir)
    snr_datat=[];
for m=1:mlen
    % NOW LOADING FROM tdatenum ONWARDS
    tdatenumt=tdatenum-mod(tdatenum,1)+m-1;
    curdtt=datetime(tdatenumt,'convertfrom','datenum');
    strdayt=char(datetime(curdtt,'format','DDD'));
    stryrst=char(datetime(curdtt,'format','yy'));
    if exist([char(snrdir(ll)),'/',num2str(tdatenumt),'.mat'],'file')==2
        load([char(snrdir(ll)),'/',num2str(tdatenumt),'.mat'])
    else
        disp('missing data')
        miss=1;
        break
    end
    snr_data(:,9)=tdatenum-mod(tdatenum,1)+m-1+snr_data(:,4)./86400;
    snr_data(:,10)=ll;
    % now get rid of data outside window
    tmpout=snr_data(:,9)<tdatenum | snr_data(:,9)>=tdatenum+tlen;
    snr_data(tmpout,:)=[];
    snr_datat=[snr_datat;snr_data];
end
if exist('miss')~=0
    break
end
snr_data_all=[snr_data_all;sortrows(snr_datat,1)];
clear snr_datat snr_data
end

if numel(snr_data_all)==0
    disp('no data - continue')
    return
end

% HERE IS WHERE THE NEW FUNCTION NEEDS TO GO
snr_dt=[];
rh_stats=[];
for ll=1:numel(snrdir)
    [rh_statst,snr_dtt] = snr2arcs(snr_data_all,staxyz,elvlims,azilims,rhlims,...
        dt,satconsts,sig,arclims,snrfigs,lspfigs,gfresnel);
    rh_statst(:,12)=ll;
    snr_dtt(:,5)=ll;
    rh_stats=[rh_stats;rh_statst];
    snr_dt=[snr_dt;snr_dtt];
end

%scatter(rh_stats(:,1),-rh_stats(:,2))
%return

if min(rh_stats(:,1))>tdatenum+tlen/3 || size(rh_stats,1)<2
    disp('not enough data - continue')
    return
end

% for adjusting heights
if numel(snrdir)>1
    for ll=1:numel(snrdir)
        in=rh_stats(:,12)==ll;
        meanhgts(ll)=nanmean(rh_stats(in,2));
        if ll>1
            rh_stats(in,2)=rh_stats(in,2)+meanhgts(1)-meanhgts(ll);
        end
    end
else
    meanhgts=0;
end

rh_stats=sortrows(rh_stats,1);
hsmooth=smoothdata(rh_stats(:,2),'movmean',5);
diff1=abs(hsmooth-rh_stats(:,2));
std1=std(diff1);
delete=diff1(:,1)>stdfac*std1; % better to keep at 3
rh_stats(delete,:)=[];

if largetides==0
indt=snr_dt(:,3)>tdatenum+tlen/3 & snr_dt(:,3)<tdatenum+2*tlen/3;
t_allt=snr_dt(indt,3);
maxtgap=max(diff(sort(t_allt)));
else 
indt=rh_stats(:,1)>tdatenum+tlen/3 & rh_stats(:,1)<tdatenum+2*tlen/3;
indt=find(indt);
if min(indt)~=1
    indt=[min(indt)-1;indt];
end
if max(indt)~=size(rh_stats,1)
    indt=[indt;max(indt)+1];
end
tempt=rh_stats(indt,1);
maxtgap=max(diff(sort(tempt)));
end
disp(['max gap is ',num2str(maxtgap*24*60),' minutes'])

if  numel(maxtgap)==0 || maxtgap>kspac*1.5 ||...
        max(rh_stats(:,1))<tdatenum+2*tlen/3 || min(rh_stats(:,1))>tdatenum+tlen/3
    disp('gap in data bigger than node spacing')
    % CHOOSE THIS
    %disp('continue with risk of instabilities')
    % OR THIS
    return
end

if size(rh_stats,1)<2
    disp('not enough data - continue')
    return
end

knots=[tdatenum*ones(1,p) ...
    tdatenum:kspac:tdatenum+tlen ...
    (tdatenum+tlen)*ones(1,p)];
nsfac=tlen/kspac+p;
sfacs_0=mean(rh_stats(:,2))*ones(1,nsfac);
tempfun_init=@(sfacs) bspline_spectral(sfacs,p,knots,rh_stats(:,4),rh_stats(:,1),1)-rh_stats(:,2).';
%tempfun_init=@(sfacs) bspline_spectral_accel(sfacs,p,knots,rh_stats(:,4),rh_stats(:,1),rh_stats(:,11),1)-rh_stats(:,2).';
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off'); % off for now but maybe should check sometimes
sfacs_init=lsqnonlin(tempfun_init,sfacs_0,[],[],options); %lsqnonlin or fsolve??
in=rh_stats(:,1)>tdatenum+tlen/3 & rh_stats(:,1)<tdatenum+2*tlen/3;
rh_stats=rh_stats(in,:);
xinit=rh_stats(:,1);
hinit=rh_stats(:,2);
xinit=xinit.';
hinit=hinit.';
sfacspre=sfacs_init;

if skipjs==0
disp('joakim part')
doprev=0;
if exist('sfacsjs')==0 || doprev==0
if largetides==1
sfacs_0=sfacs_init;
else
sfacs_0=median(rh_stats(:,2))*ones(size(sfacs_init));
end
else
inds=tlen/(3*kspac)+2;
sfacs_0=sfacsjs(inds:end);
sfacs_0=[sfacs_0(1)*ones(1,p-1) sfacs_0 sfacs_0(end)*ones(1,tlen/(3*kspac)+p-1)];
end
consts=gps+glo+gal;
sfacs_0=[sfacs_0 zeros(1,consts*2)];
if numel(roughin)>0
sfacs_0=[sfacs_0 roughin];
tempfun=@(sfacs) bspline_js(sfacs,snr_dt(:,3),snr_dt(:,1),snr_dt(:,2),knots,...
    p,snr_dt(:,4),gps,glo,gal,snr_dt(:,5),meanhgts,sig);
else
tempfun=@(sfacs) bspline_jsnorough(sfacs,snr_dt(:,3),snr_dt(:,1),snr_dt(:,2),knots,...
    p,snr_dt(:,4),gps,glo,gal,snr_dt(:,5),meanhgts,sig);
end
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off');
tic
sfacs_ls=lsqnonlin(tempfun,sfacs_0,[],[],options); %lsqnonlin or fsolve??
toc
disp('****least squares done****')

if numel(roughin)>0
sfacsjs=sfacs_ls(1:end-consts*2-1);
consts_out=sfacs_ls(end-consts*2:end-1);
roughout=sfacs_ls(end);  
else
sfacsjs=sfacs_ls(1:end-consts*2);
consts_out=sfacs_ls(end-consts*2:end);
end

else
sfacsjs=sfacspre;
end


end
