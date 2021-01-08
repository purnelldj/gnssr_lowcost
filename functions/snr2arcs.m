function [rh_stats,snr_dt] = snr2arcs(snr_data,staxyz,elvlims,azilims,rhlims,dt,satconsts,sig,arclims,snrfigs,lspfigs,gfresnel)

pwdstr=pwd;
if satconsts(2)
load([pwdstr,'/functions/glonasswlen.mat'])
end

if snrfigs==1 || lspfigs==1 || gfresnel==1
width = 3.6;     % Width in inches % 3.5 was for putting two on one line i think
height = 1.3; % was 1.3    % Height in inches
fsz = 11;      % Fontsize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*70]);
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2; %
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all
addpath([pwdstr,'/functions/fresnel'])
delete([pwdstr,'/tempoutput/*'])
end

% changing the temporal resolution
modspl=mod(snr_data(:,4),dt);
delt=modspl(:)>0;
snr_data(delt,:)=[];

% gettting rid of data outside azi and elevation lims
in=snr_data(:,3)>azilims(1) & snr_data(:,3)<azilims(2);
snr_data=snr_data(in,:);
in=snr_data(:,2)<elvlims(2) & snr_data(:,2)>elvlims(1);
snr_data=snr_data(in,:);
%if exist('azi_mask')==1
%    out=snr_data(:,3)>azi_mask(1) & snr_data(:,3)<azi_mask(2);
%    snr_data(out,:)=[];
%end
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

ind=1;
cursat=snr_data(1,1);
stind=1;
stopp=1;
rhind=0;
snr_dt=zeros(0,4);
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
    if ind-stind>1
    if snr_data(ind,1)~=cursat || ind==size(snr_data,1) || abs(curdt)>=3*dt/86400 ...
            || fwd2~=fwd1 || snr_data(ind,2)==snr_data(ind-1,2) || ind-stopp-stind>arclims(2)/dt
        if ind-stopp-stind<arclims(1)/dt
            snr_data(stind:ind-stopp,:)=[];
            ind=stind;
        else
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
        if numel(sinelvt)>2
        sinelvt=sind(sinelvt);
        snrtmp(del,:)=[];
        ttmp(del,:)=[];
        p1=polyfit(sinelvt,snrtmp,2);
        y1=polyval(p1,sinelvt);
        snrtmp=snrtmp-y1;
        snrtmp=snrtmp*100/(max(abs(snrtmp))); % NORMALIZE
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
        maxf1=numel(sinelvt)/(2*(max(sinelvt)-min(sinelvt)));
        ovs=round(Lcar/(2*prec*(max(sinelvt)-min(sinelvt))));
        if sum(diff(sinelvt)>0)==numel(sinelvt)-1 || sum(diff(sinelvt)>0)==0
        [psd,f]=plomb(snrtmp,sinelvt,maxf1,ovs,'normalized'); %
        reflh1=f.*0.5*Lcar;
        [~,id]=max(psd(:));
        if reflh1(id) > rhlims(1) && reflh1(id) < rhlims(2)
        snr_datatmp=snr_data(stind:ind-stopp,:);    
        rhind=rhind+1;
        rh_stats(rhind,1)=mean(ttmp);                            % datenum
        rh_stats(rhind,2)=reflh1(id);                            % rh (m)
        rh_stats(rhind,3)=satnotmp;                              % sat prn
        rh_stats(rhind,4)=tand(mean(snr_datatmp(:,2)))/(((pi/180)*...
            (snr_datatmp(end,2)-snr_datatmp(1,2)))...
            /((snr_datatmp(end,9)-snr_datatmp(1,9))*86400));     % tan(th)/dth/dt
        rh_stats(rhind,5)=min(snr_datatmp(:,2));                 % THETA MIN
        rh_stats(rhind,6)=max(snr_datatmp(:,2));                 % THETA MAX
        rh_stats(rhind,7)=nanmean(snr_datatmp(:,3));             % MEAN AZI
        rh_stats(rhind,8)=nanmean(y1);                           % mean mag. tSNR
        rh_stats(rhind,9)=max(psd);                              % the peak
        rh_stats(rhind,10)=var(snrtmp);                          % the variance
        rh_stats(rhind,11)=snr_datatmp(end,9)-snr_datatmp(1,9);  % arc length s
        if rh_stats(rhind,3)<0
            rh_stats(rhind,11)=-rh_stats(rhind,11);
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        % SNR PLOT
        if snrfigs==1
        figure('visible','off')
        plot(asind(sinelvt),snrtmp,'b','linewidth',0.8)
        xlabel('Elevation angle','interpreter','latex','fontsize',fsz)
        ylabel('$\delta SNR$','interpreter','latex','fontsize',fsz)
        set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
        axis([elvlims(1) elvlims(2) -120 120])
        print([pwdstr,'/tempoutput/SNR_',num2str(satnotmp),'_',...
            num2str(round(mean(ttmp))),...
            '.png'],'-dpng','-r300')
        end
        % PERIODOGRAM
        if lspfigs==1
        figure('visible','off')
        plot(reflh1,psd,'r')
        axis([0 30 0 30])
        xlabel('Reflector height','interpreter','latex','fontsize',fsz)
        ylabel('Power','interpreter','latex','fontsize',fsz)
        set(gca,'ticklabelinterpreter','latex','fontsize',fsz-1)
        print([pwdstr,'/tempoutput/LSP_',num2str(satnotmp),'_',...
            num2str(round(mean(ttmp))),...
            '.png'],'-dpng','-r300')
        end
        % FRESNEL ZONES
        if gfresnel==1
            % currently for L1 signal, could change to L2, see function description
            fresout=[pwdstr,'/tempoutput'];
            googlefresnel(lat,lon,asind(sinelvt(1)),snr_datatmp(1,3),satnotmp,fresout,reflh1(id),1)
            googlefresnel(lat,lon,asind(sinelvt(end)),snr_datatmp(end,3),satnotmp,fresout,reflh1(id),1)
        end
        %%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        end
        end
        end
        % now onto next bit of data
        end
        stind=ind;
        if ind<size(snr_data,1)
            cursat=snr_data(ind,1);
        end
    end
    end
    fwd2=fwd1;
end

end

