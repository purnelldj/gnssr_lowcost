function [sfacsjs,sfacspre,rh_stats,consts_out,roughout] = invsnr(...
    tdatenum,snrdir,staxyz,elvlims,azilims,rhlims,kspac,tlen,dt,satconsts,...
    sig,arclim,pktnlim,largetides,roughin,skipjs, meanhgt_input, smoothqc)


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
% arclim: time limit (in seconds) for spectral analysis of satellite arc
% pktnlim: for outlier detection - max power of reflector height must be
% larger than pktnlim times the mean of the periodogram outside rhlims
% largetides: changes the initial guess for node values (1 or 0)
% roughin: initial sfc roughness guess ('s') in metres for least squares
% adjustment, set to '' if don't want to use roughness
% skipjs: skip the least squares adjustment as per strandberg et. al and
% instead just fit the b-spline curve using spectral analysis estimates

% OUTPUTS
% these outputs are to go with the function 'invsnr_plot.m'
% sfacsjs: node values estimated from Strandberg et al. analysis
% sfacspre: node values estimated by using spectral analysis adjustment
% rh_stats: see output from snr2arcs
% consts_out: the other variables estimated as part of least squares
% adjustment
% roughout: the roughness parameter that is estimated in least squares
% adjustment

p=2; % bspline order

snrfigs=0;
lspfigs=0;
gfresnel=0;

gps=satconsts(1);
glo=satconsts(2);
gal=satconsts(3);

sfacsjs=NaN;
sfacspre=NaN;
rh_stats=NaN;
consts_out=NaN;
roughout=NaN;

pwdstr=pwd;
addpath([pwdstr,'/matlab_functions/'])
addpath([pwdstr,'/matlab_functions/bspline'])

tdatenum=tdatenum-tlen/3;
 
curdt=datetime(tdatenum+tlen/3,'convertfrom','datenum');
disp(char(curdt))

% work out if need two days or just one
mlen=1;
%if tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))>0
%    mlen=tdatenum+tlen-mod(tdatenum+tlen,1)-(tdatenum-mod(tdatenum,1))+1;
%    if mod(tdatenum+tlen,1)==0
%        mlen=mlen-1;
%    end
%end
dt_s = datetime(tdatenum,'convertfrom','datenum');
dt_e = datetime(tdatenum+tlen,'convertfrom','datenum');
if day(dt_s) ~= day(dt_e)
    %mlen = days(dt_e - dt_s) - mod(days(dt_e - dt_s),1) + 1;
    mlen = (tdatenum+tlen)-mod(tdatenum+tlen,1) - (tdatenum - mod(tdatenum,1)) + 1;
    if hour(dt_e) == 0
        mlen = 1;
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
normalize=1;
snr_dt=[];
rh_stats=[];
for ll=1:numel(snrdir)
    tfilter = snr_data_all(:, 10) == ll;
    snrdatat = snr_data_all(tfilter, :);
    [rh_statst,snr_dtt] = snr2arcs(snrdatat,staxyz,elvlims,azilims,...
        rhlims,dt,satconsts,sig,arclim,normalize,pktnlim,snrfigs,lspfigs,...
        gfresnel);
    rh_statst(:,12)=ll;
    snr_dtt(:,5)=ll;
    rh_stats=[rh_stats;rh_statst];
    snr_dt=[snr_dt;snr_dtt];
end

%scatter(rh_stats(:,1),-rh_stats(:,2))
%return

if min(rh_stats(:,1))>tdatenum+tlen/3 || max(rh_stats(:,1))<tdatenum+2*tlen/3 ...
        || size(rh_stats,1)<2
    disp('not enough data - continue')
    return
end

% for adjusting heights
if numel(snrdir)>1
    for ll=1:numel(snrdir)
        in=rh_stats(:,12)==ll;
        if meanhgt_input
            meanhgts = meanhgt_input;
        else
            meanhgts(ll)=nanmean(rh_stats(in,2));
        end
        if ll>1
            rh_stats(in,2)=rh_stats(in,2)+meanhgts(1)-meanhgts(ll);
        end
    end
else
    meanhgts=0;
end

std_consec=std(diff(rh_stats(:,2)));
% disp(['standard deviation is ',num2str(std_consec)])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OUTLIER DETECTION OPTIONS

%if std_consec > 3
%    disp('std bigger than chosen value')
%    return    
%end

% get rid of satellite ascending and descending points separately
if smoothqc
%presm = size(rh_stats, 1);
%disp([num2str(presm), ' points pre smooth'])
fwds=rh_stats(:,4)>0;
fwd_ids=find(fwds);
bkd_ids=find(~fwds);
rhsmooth_fwd=smoothdata(rh_stats(fwds,2),'movmean',5);
rhsmooth_bkd=smoothdata(rh_stats(~fwds,2),'movmean',5);
diffs_fwd=abs(rhsmooth_fwd-rh_stats(fwds,2));
diffs_bkd=abs(rhsmooth_bkd-rh_stats(~fwds,2));
del1=diffs_fwd>3*std(diffs_fwd);
del2=diffs_bkd>3*std(diffs_bkd);
%disp(['standard deviations are: ', num2str(std(diffs_fwd)),...
%    ' and ', num2str(std(diffs_bkd))])
fwd_ids(del1)=[];
bkd_ids(del2)=[];
all_ids=[fwd_ids;bkd_ids];
rh_stats=rh_stats(sort(all_ids),:);
%disp(['smooth got rid of ', num2str(presm - size(rh_stats,1)), ' points'])
end

% for testing dbscan algorithm to remove outliers
dbscantest=0;
if dbscantest==1
epsilon = 1.5; % search radius - adjust for tidal range
minpoints = 3;
disp('dbscale is for rotg')
%dbscale = 15; % VLIS
dbscale = 22.55; % ROTG
rht_scaled = rh_stats(:,1).*dbscale;
rh_t = rh_stats(:,2);
idx = dbscan([rht_scaled rh_t],epsilon,minpoints);
del1 = idx==-1;
rh_stats(del1,:)=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

if  numel(maxtgap)==0 || maxtgap>kspac %||...
        %max(rh_stats(:,1))<tdatenum+2*tlen/3 || min(rh_stats(:,1))>tdatenum+tlen/3
    disp('gap in data bigger than node spacing')
    % CHOOSE THIS
    disp('continue with risk of instabilities')
    % OR THIS
    %return
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
tempfun_init=@(sfacs) bspline_spectral(sfacs,p,knots,rh_stats(:,4),...
    rh_stats(:,1),1)-rh_stats(:,2).';
%tempfun_init=@(sfacs) bspline_spectral_accel(sfacs,p,knots,rh_stats(:,4),...
%rh_stats(:,1),rh_stats(:,11),1)-rh_stats(:,2).';
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off'); % off for now but maybe should check sometimes
sfacs_init=lsqnonlin(tempfun_init,sfacs_0,[],[],options); %lsqnonlin or fsolve??
in=rh_stats(:,1)>tdatenum+tlen/3 & rh_stats(:,1)<tdatenum+2*tlen/3;
rh_stats=rh_stats(in,:);
sfacspre=sfacs_init;
% disp([num2str(numel(hinit)),' points'])

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

