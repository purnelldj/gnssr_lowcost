clear all; close all;

%choose stard and end points
styr=15; % 23 jun '09 to 04 apr '10
endyr=15;
stmn=2;
endmn=2;
sday=14; % sep 18 ish for riv3
eday=17; % up to last day
station='pbay';
belly=0;
mac=1;
nodefrq=6/24; % CHANGE HOUR/24
bspline_order=2;

%datevecs=datenum(2000+styr,stmn,sday,0,0,0);
%datevece=datenum(2000+endyr,endmn,eday,0,0,0);
datevecs=0;
datevece=3;
n_control=(datevece-datevecs)/nodefrq-1;
knots=[datevecs*ones(1,bspline_order) ...
    linspace(datevecs,datevece,n_control-bspline_order+1) ...
    datevece*ones(1,bspline_order)]
return

if belly==0
    name='dave';
    if mac==0
    name='davidpurnell';
    end
datastr=strcat('/Users/',name,'/data/');
workstr=strcat('/Users/',name,'/Library/Mobile Documents/com~apple~CloudDocs/work/gnssr_matlab/');
else
datastr='/home/dpurnell/projects/ctb-ng50/dpurnell/data/';
workstr='/home/dpurnell/projects/ctb-ng50/dpurnell/work/gnssr_matlab/';
end

if strcmp(station(1:4),'sc02')==1
    load(strcat(datastr,'/sc02/tg_2015_6min.mat'))
elseif strcmp(station,'pbay')==1
    load(strcat(datastr,station,'/tg_2015_6min.mat'))
elseif strcmp(station(1:4),'sab2')==1
    load(strcat(datastr,'sab1/tg_feb19_mar10.mat'))
end

plotl=6/(24*60);
xinvt=datevecs:plotl:datevece;
slvl2=interp1(xaxis,slvl,xinvt,'nearest');
% FIRST NEED TO ESTIMATE CONTROL POINTS BASED ON INITIAL H GUESS
coefs_0=mean(slvl2)*ones(1,n_control);
tempfun=@(coefs) bspline_deboor(bspline_order+1,knots,coefs,xinvt)-slvl2;
options=optimoptions(@lsqnonlin,'Algorithm','levenberg-marquardt',...
    'Display','off');
coefs_ls=lsqnonlin(tempfun,coefs_0,[],[],options); %lsqnonlin or fsolve??
slvlspline=bspline_deboor(bspline_order+1,knots,coefs_ls,xinvt);

plot(xinvt,slvl2,'k')
hold on
plot(xinvt,slvlspline,'r')

rmss=rms(slvl2-slvlspline);
disp(['RMS is ',num2str(rmss*100),' cm'])



