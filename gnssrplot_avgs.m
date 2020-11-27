clear all
close all

% to plot GNSS_R sea level measurements as per Larson et al. (2017)
% peaks from LSPs are CORRECTED for the dynamic sea surface by fitting
% tidal constituents to the data via least squares
% the inputs are:
% '.mat' files with the variable 'slvlr'
% 'slvlr' contains information from real or fake GNSS-R data
% the outputs are:
% rms_adj = RMS of CORRECTED GNSS-R sea level measurements
% rms_nadj = RMS of GNSS-R sea level measurements prior to correction
% rh_adj = time series of CORRECTED GNSS-R sea level measurements
% rh_nadj = time series of GNSS-R sea level measurements prior to correction
% t_rh = timings of GNSS-R sea level measurements
% to run this code you need:
% a 'station_inputs' directory
% data contained in 'datastr' to plot
% tide gauge data or model tides
datenums=datenum(2020,5,17,0,0,0);
datenume=datenum(2020,6,1,0,0,0);

% change variables
stations=cellstr({'sab2a'});
plotnow=1; % choose whetherr or not to generate plot
l1=1; % L1 signal choose (1) or (0)
l2=0; % L2 signal choose (1) or (0)
gps=1; % choose (1) or (0)
glo=1; % choose (1) or (0)
gal=1; % choose (1) or (0)
% for looking at the direction of satellite overpassess
justfwd=0; % 1 for fwd, 2 for not fwd, 0 for off (all data)
mac=1; % change local directory
belly=0;
doelvlims=0;
autstr='';
remove_mean=0;
tidegauge=1;
histogramm=0;
hourlyavg=1;
remvoutliers=1;


% data drectory
dirstr='obs_stats_1050';

% add your local directory to datastr
if mac==1
    name='dave';
    datastr=strcat('/Users/',name,'/data/');
elseif mac==0
    name='davidpurnell';
    datastr=strcat('/Users/',name,'/data/');
end
if belly==1
    datastr='/home/dpurnell/projects/ctb-ng50/dpurnell/data/';
    %workstr='/home/dpurnell/projects/ctb-ng50/dpurnell/work/gnssr_matlab/';
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if datenums>=datenume
    disp('you have not chosen your dates well')
    return
end

% HERE IS WHERE YOU ARE LOADING ALL DATA TO BE ANALYZED
slvlr_all=[];
tdatenum=datenums-1;
while tdatenum<datenume
    tdatenum=tdatenum+1;
start=1;
for bb=1:numel(stations)
    station=char(stations(bb));
if exist(strcat(datastr,station,'/',dirstr,'/',num2str(tdatenum),'.mat'),'file')==2
    load(strcat(datastr,station,'/',dirstr,'/',num2str(tdatenum),'.mat'))
    slvlr(:,15)=bb;
    slvlr_all=[slvlr_all;slvlr];
end
end
end
slvlr=slvlr_all;

if size(slvlr,1)<1
    disp('apparently there are no data points in this time range')
    return
end

% delete outside of time frame
out=slvlr(:,1)<datenums | slvlr(:,1)>datenume;
slvlr(out,:)=[];
% fwd delete
if justfwd==1
out=slvlr(:,3)<0;
slvlr(out,:)=[];
elseif justfwd==2
out=slvlr(:,3)>0;
slvlr(out,:)=[];
end

% getting desired satellite constellations
if glo==0
    delete=slvlr(:,2)>32 & slvlr(:,2)<32+24;
    slvlr(delete,:)=[];
end
if gps==0
    delete=slvlr(:,2)<33;
    slvlr(delete,:)=[];
end
if gal==0
    delete=slvlr(:,2)>32+24;
    slvlr(delete,:)=[];
end

% getting rid of bad elv lims
if doelvlims==1
    %run(strcat('station_inputs/',station,'_input',autstr))
    %clear makesnr
    elv_low=10;
    elv_high=20;
    elvlims=1;
    delete=abs(slvlr(:,4)-elv_low)>elvlims;
    slvlr(delete,:)=[];
    delete=abs(slvlr(:,5)-elv_high)>elvlims;
    slvlr(delete,:)=[];
end

% GET RID OF NANS
if l1==1 && l2==1
keep=~isnan(slvlr(:,7));
tsize=size(slvlr,1)+1;
slvlr=[slvlr;slvlr(keep,:)];
slvlr(tsize:end,7)=slvlr(tsize:end,11);
slvlr(:,11)=NaN;
elseif l1==0 && l2==1
keep=isnan(slvlr(:,11))==0;
slvlr=slvlr(keep,:);
slvlr(:,7)=slvlr(:,11);
end
delete=isnan(slvlr(:,7))==1;
slvlr(delete,:)=[];
slvlr=sortrows(slvlr,1);

if remove_mean==1
% removing mean from time series
tempmean=nanmean(slvlr(:,7));
slvlr(:,7)=tempmean-slvlr(:,7);
tempmean=nanmean(slvlr(:,11));
slvlr(:,11)=tempmean-slvlr(:,11);
end

rh_nadj=slvlr(:,7);
t=slvlr(:,1);
statinds=slvlr(:,15);

disp(['mean height is ',num2str(nanmean(rh_nadj))])

if remvoutliers==1
rh_nadj2=[];
t2=[];
statinds2=[];
for bb=1:numel(stations)
    indt=statinds(:)==bb;
    rht=rh_nadj(indt);
    tt=t(indt);
    statindst=statinds(indt);
hsmooth=smoothdata(rht,'movmean',5);
diff1=abs(rht-hsmooth);
std1=std(diff1);
delete=diff1(:,1)>3*std1; % here choose sigma bounds to remove
rht(delete)=[];
tt(delete)=[];
statindst(delete)=[];
rh_nadj2=[rh_nadj2;rht];
t2=[t2;tt];
statinds2=[statinds2;statindst];
end
rh_nadj=rh_nadj2;
t=t2;
statinds=statinds2;
end

if hourlyavg==1
    ddtt=1/24;
    thourly=datenums+0.5*ddtt:ddtt:datenume+1-0.5*ddtt;
    rhourly=NaN(numel(stations),numel(thourly));
    for bb=1:numel(stations)
        indt=statinds(:)==bb;
        rht=rh_nadj(indt);
        tt=t(indt);
        for cc=1:numel(thourly)
            indtt=tt(:)>=thourly(cc)-0.5*ddtt & tt(:)<thourly(cc)+0.5*ddtt;
            if sum(indtt)>2
                %rhourly(bb,cc)=nanmedian(rht(indtt));
                rhourly(bb,cc)=nanmean(rht(indtt));
            else
                rhourly(bb,cc)=NaN;
            end
        end
        %checkmeans(bb)=nanmean(rhourly(bb,:));
        rhourly(bb,:)=nanmean(rhourly(bb,:))-rhourly(bb,:);
        %plot(thourly,rhourly(bb,:))
        hold on
    end
    allhourly=nanmedian(rhourly,1);
    %return
end

if tidegauge==1
    tmpsta=char(stations(bb));
    if strcmp(tmpsta(1:4),'sab1')
        load('~/data/sab1/tg_dec1_feb17.mat')
        out=xaxis(:)<datenums | xaxis(:)>datenume;
        slvl(out)=[];
        xaxis(out)=[];
    end
    if strcmp(tmpsta(1:4),'sab2')
        load('~/data/sab2/tg_jan1_jun15.mat')
        out=xaxis(:)<datenums | xaxis(:)>datenume;
        slvl(out)=[];
        xaxis(out)=[];
        tomakepaperplot=0;
        if tomakepaperplot==1
            slvl=-slvl+25.2;
        end
    end
    if strcmp(char(stations),'osou')
        load('~/data/osou/osoutg.mat')
        out=xaxis(:)<datenums | xaxis(:)>datenume;
        slvl(out)=[];
        xaxis(out)=[];
    end
    if strcmp(tmpsta(1:4),'rv3s')
        load('~/data/rv3s/tg_sep9_oct10.mat')
        out=xaxis(:)<datenums | xaxis(:)>datenume;
        slvl(out)=[];
        xaxis(out)=[];
    end
end

if hourlyavg==1
    rhnnans=~isnan(allhourly(:));
    thourly=thourly(rhnnans);
    allhourly=allhourly(rhnnans);
    tghourly=interp1(xaxis,slvl,thourly,'nearest');
    tghourly=tghourly-nanmean(tghourly);
    rmshourly=rms(tghourly-allhourly)
    plot(thourly,tghourly,'k')
    hold on
    plot(thourly,allhourly,'r')
    return
end

%return
osou_levels=0;
if osou_levels==1
    rh_nadj=-rh_nadj+(3.9562+0.15508);
    if tidegauge==1
    slvl=slvl+2.56467;
    end
end

getdiffs=0;
if getdiffs==1
    slvltg_rms=interp1(xaxis,slvl,t,'nearest');
    meansep=nanmean(rh_nadj-slvltg_rms);
    disp(['GNSS-R is ',num2str(meansep*100),' cm above the tide gauge'])
end

if plotnow==1
% PLOTTING
% Defaults
width =7;     % Width in inches % 3.5 was for putting two on one line i think
height =4;    % Height in inches
alw = 0.75;    % AxesLineWidth
fs =15;      % Fontsize
lw = 1;      % LineWidth
msz = 4;       % MarkerSize
% The properties we've been using in the figures
set(0,'defaultLineLineWidth',lw);   % set the default line width to lw
set(0,'defaultLineMarkerSize',msz); % set the default line marker size to msz
% Set the default Size for display
defpos = get(0,'defaultFigurePosition');
set(0,'defaultFigurePosition', [defpos(1) defpos(2) width*100, height*70]);
% Set the defaults for saving/printing to a file
set(0,'defaultFigureInvertHardcopy','on'); % This is the default anyway
set(0,'defaultFigurePaperUnits','inches'); % This is the default anyway
defsize = get(gcf, 'PaperSize');
left = (defsize(1)- width)/2; %
%left=0;
bottom = (defsize(2)- height)/2;
defsize = [left, bottom, width, height];
set(0, 'defaultFigurePaperPosition', defsize);
close all  

figure('visible','on')
for ij=1:numel(stations)
    indt=statinds(:)==ij;
    pointsperday=sum(indt)/(datenume-datenums);
    disp(['points per day = ',num2str(pointsperday)])
    scatter(t(indt),rh_nadj(indt),'+','linewidth',lw)
    hold on
end
if tidegauge==1
    plot(xaxis,slvl,'k')
end
%h=legend([stations,cellstr({'tide gauge'})],'numcolumns',numel(stations)+1);
h=legend([cellstr({'antenna a','antenna b','antenna c','antenna d'})],'numcolumns',4);
axis([datenums datenume -inf inf])
%set(gca,'xtick',datenums:4/24:datenume)
datetick('x',1,'keeplimits','keepticks')
box on
%xlabel('Day of the year, 2017','interpreter','latex','fontsize',fs)
ylabel('Reflector Height (m)','interpreter','latex','fontsize',fs)
set(gca,'ticklabelinterpreter','latex','fontsize',fs)
set(h,'interpreter','latex','fontsize',fs,'location','north')
print('test', '-dpng', '-r300');
end
