function [rh_stats,snr_dt] = snr2arcs(snr_data,staxyz,elvlims,azilims,...
    rhlims,dt,satconsts,sig,arclim,normalize,pktnlim,snrfigs,lspfigs,gfresnel)
%%

% This code is for organising SNR data into satellite arcs and performing
% spectral analysis
% written by Dave Purnell
% https://github.com/purnelldj/gnssr_lowcost

% INPUTS:
% snr_data: output from nmea2snr.m
% staxyz: approx xyz of statiom
% elvlims: elevation angle limits
% azilims: azimuth angle limits
% rhlims: reflector height limits min/max (m)
% dt: temporal resolution of input SNR data (in seconds), if larger than
% the actual resolution then will remove data to convert
% satconsts: 1 by 3 double, 1 or 0 to include satellite constellations
% GPS, GLONASSS, GALILEO (e.g., [1 0 1] for GPS and GALILEO)
% sig: 1 for l1 or 2 for l2
% arclim: time limit (in seconds) for spectral analysis of satellite arc
% pktnlim: for outlier detection - max power of reflector height must be
% larger than pktnlim times the mean of the periodogram outside rhlims
% snrfigs: 1 or 0 to produce figures of SNR data in tempoutput directory
% lspfigs: as above but with LSP figures
% gfresnel: output fresnel zones of start / end points of satellite arcs

% OUTPUTS
% rh_stats: organised reflector height estimates and other stats from
% satellite arcs. See below for description of data
% snr_dt: detrended SNR data for inverse analysis (e.g., invsnr.m)

rh_stats=[];
snr_dt=[];

pwdstr=pwd;
if satconsts(2)
load([pwdstr,'/matlab_functions/glonasswlen.mat'])
end

addpath([pwdstr,'/matlab_functions/fresnel'])

if snrfigs==1 || lspfigs==1 || gfresnel==1
width = 4.5;     % Width in inches % 3.5 was for putting two on one line i think
height = 2; % was 1.3    % Height in inches
fsz = 11;      % Fontsize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*70]);
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2; %
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all
delete([pwdstr,'/tempoutput/*'])
end

% GET RID OF SAT E36 - doesn't exist??
tdel = snr_data(:, 1) == 92;
snr_data(tdel, :) = [];

% changing the temporal resolution
modspl=mod(snr_data(:,4),dt);
delt=modspl(:)>0;
snr_data(delt,:)=[];

% gettting rid of data outside azi and elevation lims
in=snr_data(:,3)>azilims(1) & snr_data(:,3)<azilims(2);
snr_data=snr_data(in,:);
in=snr_data(:,2)<elvlims(2) & snr_data(:,2)>elvlims(1);
snr_data=snr_data(in,:);

% change units from db-hz to watt/watt
tmp=snr_data(:,7)==0;
snr_data(tmp,7)=NaN;
tmp=snr_data(:,8)==0;
snr_data(tmp,8)=NaN;
snr_data(:,7)=sqrt(10.^(snr_data(:,7)./10));
snr_data(:,8)=sqrt(10.^(snr_data(:,8)./10));

% getting desired satellite constellations
if ~satconsts(1)
    del=snr_data(:,1)<32+1;
    snr_data(del,:)=[];
end
if ~satconsts(2)
    del=snr_data(:,1)>32 & snr_data(:,1)<32+24+1;
    snr_data(del,:)=[];
end
if ~satconsts(3)
    del=snr_data(:,1)>32+24;
    snr_data(del,:)=[];
end

if gfresnel==1
    [lat,lon,~]=ecef2lla(staxyz(1),staxyz(2),staxyz(3));
    lat=lat/pi*180;
    lon=lon/pi*180;
    lon(lon>180)=lon-360;
end

snr_data=sortrows(snr_data,1);

ind=1;
cursat=snr_data(1,1);
stind=1;
stopp=1;
rhind=0;
snr_dt=zeros(0,4);
%rh_stats=zeros(0,11);
prec=0.001; % precision of lomb scargle
if snr_data(2,2)-snr_data(1,2)<0
    fwd2=0;
else
    fwd2=1;
end
while stopp==1
    ind=ind+1;
    if ind==size(snr_data,1)
        stopp=0;
    end
    if snr_data(ind,2)-snr_data(ind-1,2)<0
        fwd1=0;
    else
        fwd1=1;
    end
    curdt=snr_data(ind,9)-snr_data(ind-1,9);
    curdtarc = snr_data(ind-stopp, 9) - snr_data(stind, 9);
    if ind-stind>1
    if snr_data(ind,1)~=cursat || ind==size(snr_data,1) || abs(curdt) > 1/(24*12) ...
            || fwd2~=fwd1 || snr_data(ind,2)==snr_data(ind-1,2) || curdtarc > (arclim-dt/2)/86400 
        sinelvt=snr_data(stind:ind-stopp,2);
        if sig==1
        snrtmp=snr_data(stind:ind-stopp,7);
        elseif sig==2
        snrtmp=snr_data(stind:ind-stopp,8);
        end
        ttmp=snr_data(stind:ind-stopp,9);
        satnotmp=snr_data(stind,1);
        %%%%%%%%%%%%%%%%%%%
        del=isnan(snrtmp(:,1));
        sinelvt(del,:)=[];
        if numel(sinelvt)>20 % 10 is chosen arbitrarily to avoid bad fitting
        sinelvt=sind(sinelvt);
        snrtmp(del,:)=[];
        ttmp(del,:)=[];
        p1=polyfit(sinelvt,snrtmp,2);
        y1=polyval(p1,sinelvt);
        snrtmp=snrtmp-y1;
        if normalize==1
        snrtmp=snrtmp*100/(max(abs(snrtmp))); % NORMALIZE
        end
        %%%%%%%
        %sinelv=[sinelv;sinelvt];
        %snr=[snr;snrtmp];
        %t=[t;ttmp];
        %satno=[satno;satnotmp*ones(size(ttmp))];
        snr_dt=[snr_dt; sinelvt snrtmp ttmp satnotmp*ones(size(ttmp))];
        % TO GET INITIAL H TIME SERIES
        if sinelvt(2,1)-sinelvt(1,1)<0
            sinelvt=flipud(sinelvt);
            snrtmp=flipud(snrtmp);
        end
        if sig==1
        if snr_data(stind,1)<33
        Lcar=299792458/1575.42e06;
        elseif snr_data(stind,1)>32 && snr_data(stind,1)<57
        Lcar=glonassl1(snr_data(stind,1)-32);
        elseif snr_data(stind,1)>56
        Lcar=299792458/1575.42e06;
        end
        elseif sig==2
        if snr_data(stind,1)<33
        Lcar=299792458/1227.60e06;
        elseif snr_data(stind,1)>32 && snr_data(stind,1)<57
        Lcar=glonassl2(snr_data(stind,1)-32);
        elseif snr_data(stind,1)>56
        Lcar=299792458/1227.60e06;
        end 
        end
        if (sum(diff(sinelvt)>0)==numel(sinelvt)-1 || sum(diff(sinelvt)>0)==0) &&...
                 (curdtarc > (arclim-dt/2)/86400 || arclim == inf) 
            % first condition above is to make sure that points are all
            % going the same direction
        %maxf1=2*(rhlims(2)+rhlims(1))/Lcar;
        %ovs=round(Lcar/(2*prec*(max(sinelvt)-min(sinelvt))));
        %[psd,f]=plomb(snrtmp,sinelvt,maxf1,ovs,'normalized'); %
        maxf = 2 * (rhlims(1) + rhlims(2)) / Lcar;
        precisionf = 2 * 0.001 / Lcar;
        f_in = precisionf:precisionf:maxf;
        [psd,f]=plomb(snrtmp,sinelvt,f_in); %
        reflh1=f.*0.5*Lcar;
        inbounds=reflh1(:)>rhlims(1) & reflh1(:)<rhlims(2);
        psdt=psd(inbounds);
        reflht=reflh1(inbounds);
        ft = f(inbounds);
        [~,id]=max(psdt(:));
        pktn=psdt(id)/mean(psd(~inbounds));
        %if reflh1(id) > rhlims(1) && reflh1(id) < rhlims(2)
        if pktn > pktnlim && id ~=0 && id ~= numel(ft)
        snr_datatmp=snr_data(stind:ind-stopp,:);    
        rhind=rhind+1;
        rh_stats(rhind,1)=mean(ttmp);                            % datenum
        rh_stats(rhind,2)=reflht(id);                            % rh (m)
        rh_stats(rhind,3)=satnotmp;                              % sat prn
        rh_stats(rhind,4)=tand(mean(snr_datatmp(:,2)))/(((pi/180)*...
            (snr_datatmp(end,2)-snr_datatmp(1,2)))...
            /((snr_datatmp(end,9)-snr_datatmp(1,9))*86400));     % tan(th)/dth/dt
        rh_stats(rhind,5)=min(snr_datatmp(:,2));                 % THETA MIN
        rh_stats(rhind,6)=max(snr_datatmp(:,2));                 % THETA MAX
        rh_stats(rhind,7)=nanmean(snr_datatmp(:,3));             % MEAN AZI
        rh_stats(rhind,8)=nanmean(y1);                           % mean mag. tSNR
        rh_stats(rhind,9)=max(psdt);                             % the peak
        rh_stats(rhind,10)=var(snrtmp);                          % the variance
        rh_stats(rhind,11)=(snr_datatmp(end,9)-snr_datatmp(1,9))*86400;  % arc length s
        if rh_stats(rhind,3)<0
            rh_stats(rhind,11)=-rh_stats(rhind,11);
        end
        rh_stats(rhind,12)=pktn;                                 % peak to noise
        rh_stats(rhind,13)=size(snr_datatmp,1);
        rh_stats(rhind,14)=numel(f);
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        % SNR PLOT
        if snrfigs==1
        figure('visible','off')
        plot(asind(sinelvt),snrtmp,'b','linewidth',0.8)
        xlabel('Elevation angle','interpreter','latex','fontsize',fsz)
        ylabel('$\delta SNR$ (Watt/Watt)','interpreter','latex','fontsize',fsz)
        set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
        axis([elvlims(1) elvlims(2) -220 220])
        print([pwdstr,'/tempoutput/SNR_',num2str(satnotmp),'_',...
            num2str(hour(datetime(mean(ttmp), 'convertfrom', 'datenum'))),'_'...
            num2str(minute(datetime(mean(ttmp), 'convertfrom', 'datenum'))),...
            '.png'],'-dpng','-r300')
        end
        % PERIODOGRAM
        if lspfigs==1
        figure('visible','off')
        plot(reflh1,psd,'r')
        axis([0 rhlims(1)+rhlims(2) 0 inf])
        xlabel('Reflector height','interpreter','latex','fontsize',fsz)
        ylabel('Power','interpreter','latex','fontsize',fsz)
        set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
        print([pwdstr,'/tempoutput/LSP_',num2str(satnotmp),'_',...
            num2str(hour(datetime(mean(ttmp), 'convertfrom', 'datenum'))),'_'...
            num2str(minute(datetime(mean(ttmp), 'convertfrom', 'datenum'))),...
            '.png'],'-dpng','-r300')
        end
        % FRESNEL ZONES
        if gfresnel==1
            % currently for L1 signal, could change to L2, see function description
            fresout=[pwdstr,'/tempoutput'];
            googlefresnel(lat,lon,asind(sinelvt(1)),snr_datatmp(1,3),satnotmp,...
                fresout,reflh1(id),1)
            googlefresnel(lat,lon,asind(sinelvt(end)),snr_datatmp(end,3),satnotmp,...
                fresout,reflh1(id),1)
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        end
        end
        end
        % now onto next bit of data
        %end
        stind=ind;
        if ind<size(snr_data,1)
            cursat=snr_data(ind,1);
        end
    end
    end
    fwd2=fwd1;
end

rh_stats=sortrows(rh_stats,1);

end

