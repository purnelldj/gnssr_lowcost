function [t_rh,rh_invjs,rh_invpre,rh_stats,rms_js,rms_pre] = invsnr_plot(startdate,enddate,invdir,kspac,tlen,...
    plotl,tgstring,rhlims,plotrh,plotvdc)

%%

% this codes takes the output from invsnr.m and plots it and compares with
% a tide gauge (if there is one)

% INPUTS
% startdate: in datenum format
% enddate: in datenum format
% invdir: directory that contains invsnr.m output to plot
% kspac: average node spacing in days (e.g., 2/24 is 2 hours)
% tlen: length of window for analysis in days, split into 3 and the middle
% period is saved (e.g., setting to 3 means that the middle day is saved)
% plotl: the frequency of output time series, in days (e.g., 1/24 is
% hourly)
% tgstring: string of path to tide gauge data (or [] if no tide gauge data)
% rhlims: min max limits for reflector height values
% plotrh: set to 1 if you want to plot / save a figure of reflector heights
% compared with tide gauge
% plotvdc: create a van de casteele diagram

% OUPUTS
% t_rh: time vector in datenum format
% rh_invjs: output of b-spline reflector heights from Strandberg et al. analysis
% rh_invpre: output of b-spline reflector heights using spectral analysis
% rh_stats: see output from snr2arcs
% rms_js: rms between tide gauge and invjs output
% rms_pre: rms between tide gauge and invdp output

pwdstr=pwd;
addpath([pwdstr,'/matlab_functions/bspline'])

meanormedian=2;
p=2;

knots=[(startdate-tlen/3)*ones(1,p) ...
    startdate-tlen/3:kspac:enddate+1+tlen/3 ...
    (enddate+1+tlen/3)*ones(1,p)];
t_rh=startdate:plotl:enddate+1;

sfacspreall=[];
sfacsjsall=[];
rh_statsall=[];
roughness_all=[];
dispmissed=0;
for jj=1:numel(invdir)
tdatenum=startdate-tlen/3;
while round(tdatenum,10,'significant')<round(enddate+1-tlen/3,10,'significant')
    tdatenum=tdatenum+tlen/3;
    inds=tlen/(3*kspac)+2; % was +2
    inde=(2*tlen)/(3*kspac)+1; % was +1
    if tdatenum==startdate
        inds=1;
    end
    if round(tdatenum,10,'significant')>=...
            round(enddate+1-tlen/3,10,'significant')
        inde=tlen/kspac+p;
    end
    if exist([char(invdir(jj)),'/',num2str(round(tdatenum,10,'significant')),'.mat'],'file')==2
    load([char(invdir(jj)),'/',num2str(tdatenum),'.mat'])
    if (numel(sfacsjs)>1 || ~isnan(sfacsjs(1))) && size(rh_stats,1)>0
    sfacspreall=[sfacspreall sfacspre(inds:inde)];
    sfacsjsall=[sfacsjsall sfacsjs(inds:inde)];
    rh_statsall = [rh_statsall; rh_stats];
    if exist('roughness','var')==1
    roughness_all=[roughness_all roughness];
    else
    roughness_all=[roughness_all NaN];
    end
    else
        if dispmissed==0
        disp('missingdata')
        dispmissed=1;
        end
    sfacspreall=[sfacspreall NaN(1,inde-inds+1)];
    sfacsjsall=[sfacsjsall NaN(1,inde-inds+1)];
    roughness_all=[roughness_all NaN];
    end
    else
        if dispmissed==0
        disp('missingdata')
        dispmissed=1;
        end
    sfacspreall=[sfacspreall NaN(1,inde-inds+1)];
    sfacsjsall=[sfacsjsall NaN(1,inde-inds+1)];
    roughness_all=[roughness_all NaN];
    end
end    
end

rh_stats = rh_statsall;

if numel(invdir)>1
    indfac=numel(sfacsjsall)/numel(invdir);
    for kk=1:numel(invdir)
        meanhgts(kk)=nanmean(sfacsjsall((kk-1)*indfac+tlen/(3*kspac)+2:kk*indfac-tlen/(3*kspac)-1));
        meanhgts2(kk)=nanmean(sfacspreall((kk-1)*indfac+tlen/(3*kspac)+2:kk*indfac-tlen/(3*kspac)-1));
        if kk>1
            sfacsjsall((kk-1)*indfac+1:kk*indfac)=...
                sfacsjsall((kk-1)*indfac+1:kk*indfac)+meanhgts(1)-meanhgts(kk);
            sfacspreall((kk-1)*indfac+1:kk*indfac)=...
                sfacspreall((kk-1)*indfac+1:kk*indfac)+meanhgts2(1)-meanhgts2(kk);
        end
    end
    for ii=1:indfac
        jj=0;
        dptmp=[];
        jstmp=[];
        while jj<numel(invdir)
            jj=jj+1;
            dptmp=[dptmp sfacspreall(ii+(jj-1)*indfac)];
            jstmp=[jstmp sfacsjsall(ii+(jj-1)*indfac)];
        end
        if meanormedian==1
        sfacsprep(ii)=nanmean(dptmp);
        sfacsjsp(ii)=nanmean(jstmp);
        else
        sfacsprep(ii)=nanmedian(dptmp);
        sfacsjsp(ii)=nanmedian(jstmp);
        end
    end
    if numel(roughness_all)>1
    roughlen=enddate-startdate;
    for ii=1:roughlen
        jj=0;
        roughtmp=[];
        while jj<numel(invdir)
            jj=jj+1;
            roughtmp=[roughtmp roughness_all(ii+(jj-1)*roughlen)];
        end
        if meanormedian==1
        roughmean(ii)=nanmean(roughtmp);
        else
        roughmean(ii)=nanmedian(roughtmp);
        end
    end
    end
else
    sfacsprep=sfacspreall;
    sfacsjsp=sfacsjsall;
    roughmean=roughness_all;
end

% extra bit with rhlims
sfacsjsp(sfacsjsp<rhlims(1) | sfacsjsp>rhlims(2))=NaN;
sfacsprep(sfacsjsp<rhlims(1) | sfacsprep>rhlims(2))=NaN;

rh_invjs=bspline_deboor(p+1,knots,sfacsjsp,t_rh);
rh_invpre=bspline_deboor(p+1,knots,sfacsprep,t_rh);

if numel(tgstring)>1
load(tgstring)
slvl=interp1(xaxis,slvl,t_rh,'linear');
tidey=slvl;
tidex=t_rh;
rh_invjs=nanmean(rh_invjs)-rh_invjs;
rh_invpre=nanmean(rh_invpre)-rh_invpre;
rh_stats(:,2)=nanmean(rh_stats(:,2))-rh_stats(:,2);
in=isnan(tidey)==0 & isnan(rh_invjs)==0;
tidey=tidey-mean(tidey(in));
rms_js=rms(rh_invjs(in)-tidey(in));
in=isnan(tidey)==0 & isnan(rh_invpre)==0;
rms_pre=rms(rh_invpre(in)-tidey(in));
corr_js=corrcoef(rh_invjs(in),tidey(in));
corr_pre=corrcoef(rh_invpre(in),tidey(in));
disp(['rms pre is ',num2str(rms_pre*100),' cm'])
disp(['rms js is ',num2str(rms_js*100),' cm'])
disp(['correlation pre is ',num2str(corr_pre(2))])
disp(['correlation js is ',num2str(corr_js(2))])
pointsperday=sum(~isnan(rh_stats(:,2)))/(enddate+1-startdate);
disp(['points per day = ',num2str(pointsperday)])
else
rms_js=NaN;
rms_pre=NaN;
end

if plotvdc==1 || plotrh==1
width = 10;     % Width in inches % 3.5 was for putting two on one line i think
height = 5; % was 1.3    % Height in inches
fsz = 11;      % Fontsize
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*70]);
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2; %
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all
end

if plotvdc==1
figure('visible','on')
in=isnan(tidey)==0 & isnan(rh_invjs)==0;
diffs_vdc=tidey(in)-rh_invjs(in);
scatter(diffs_vdc,tidey(in),1,'.','markeredgecolor',[1.00,0.07,0.65])
xtemp=zeros(1,61);
ytemp=-30:1:30;
hold on
plot(xtemp,ytemp,'k')
axis([-1 1 min(tidey(in))+0.1 max(tidey(in))+0.1])
ylabel('Water level (m)','interpreter','latex','fontsize',fsz)
xlabel('Error, TG minus GNSS-R (m)','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz)
return
end

if plotrh==1
figure('visible','on')
if numel(tgstring)>0
plot(tidex,tidey,'k','linewidth',1.5)
hold on
end
scatter(rh_stats(:,1),rh_stats(:,2),'k+')
hold on
plot(t_rh,rh_invjs,'color',[1.00,0.07,0.65],'linewidth',1.5)
hold on
plot(t_rh, rh_invpre)
axis([min(t_rh) max(t_rh) -inf inf])
ylabel('Water level (m)','interpreter','latex','fontsize',fsz)
set(gca,'ticklabelinterpreter','latex','fontsize',fsz)
datetick('x',1,'keeplimits','keepticks')
%set(gca,'xtick',datenum(2020,9,9):0.5:datenum(2020,9,11))
h=legend('Tide gauge','GNSS-R: spectral analysis','GNSS-R: inverse method');
set(h,'interpreter','latex','fontsize',fsz)
print('invfig', '-dpng', '-r300');
end


end
