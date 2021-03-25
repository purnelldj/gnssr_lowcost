function nmea2snr(nmeafile,sp3dir,sp3option,elvlims,azilims,tempres,snrdir)

%%

% This function reads NMEA 0183 format GPS data (e.g., recorded by
% low-cost hardware) and converts it into files for GNSS-R analysis
% written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost

% NOTE: INTERPOLATING SP3 FILES PROPERLY IS IMPORTANT
%
% INPUTS
%
% nmeafile: path to nmeafile
% sp3dir: path to sp3 orbit data
% sp3option: type of sp3 data (see below) (3 is recommended)
% elvlims: (1,2) double of min and max elevation angles
% azilims: (1,2) double of min and max azimuth angles
% tempres: change temporal resolution of output data (seconds)
% snrdir: the directory which output files are saved to (see below for
% format of variables)
%
% OUTPUTS
% 
% for every day of data, a file (named as datenum of date) is saved in the
% directory 'snrdir' that contains the variables: 
% 'snr_data','staxyz','lat','lon','height'
%
% 'snr_data' is the main 8 column output
% in the same format as output from RinexSNRv2
% column 1 is satellite PRN
% column 2 is elevation angle (degrees)
% column 3 is azimuth angle (degrees)
% column 4 is time in seconds of day
% column 7 is the SNR data in dB-Hz (main output)
% 
% 'staxyz' is the mean ECEF XYZ coordinates of station (m)
% 'lat' is the mean latitude (degreees)
% 'lon' is the mean longitude (degrees)
% 'height' is the height above WGS84 ellipsoid (m)

warning('off','all')

pwdstr=pwd;
addpath([pwdstr,'/functions'])
%addpath([pwdstr,'/functions/geodetic299/geodetic'])    

donefile=0;
datevecall=[];
cnt=0;
    
while donefile==0

fid=fopen(nmeafile);
if cnt>0
    cntt=0;
    while cntt<cnt
    cntt=cntt+1;
    tline=fgets(fid);
    end
end

lat=[];
lon=[];
height=[];
s1data=NaN(86400,92,5);
disp('sorting through nmea obs')
while ~feof(fid) % this isn't working, sort out
    cnt=cnt+1;
    tline=fgets(fid);
    if size(tline,2)<5
        continue
    end
    if strcmp(tline(4:6),'RMC')==1
        tmp=textscan(tline,'%s','Delimiter',',');
        datetmp=tmp{1,1}{10,1};
        dd=str2double(datetmp(1:2));
        mm=str2double(datetmp(3:4));
        yyyy=2000+str2double(datetmp(5:6));
        if ismember(datenum(yyyy,mm,dd),datevecall)==0
            datevecall=[datevecall;datenum(yyyy,mm,dd)];
            if numel(datevecall)>1
                disp(char(datetime(datenum(yyyy,mm,dd)-1,'convertfrom','datenum')))
                break
            end
        end
        datevecn=datenum(yyyy,mm,dd);
    elseif strcmp(tline(4:6),'GGA')==1
        tmp=textscan(tline,'%s','Delimiter',',');
        ggalat=tmp{1,1}{3,1};
        if numel(ggalat)>1 && size(tmp{1,1},1)>11
        latt=str2double(ggalat(1:2)) + (str2double(ggalat(3:end))/60);
        ggalats=tmp{1,1}{4,1};
        if strcmp(ggalats,'S')==1
            latt=-latt;
        end
        lat=[lat;latt];
        ggalon=tmp{1,1}{5,1};
        lont=str2double(ggalon(1:3)) + (str2double(ggalon(4:end))/60);
        ggalonw=tmp{1,1}{6,1};
        if strcmp(ggalonw,'W')==1
            lont=-lont;
        end
        lon=[lon;lont];
        height=[height;str2double(tmp{1,1}{10,1})+str2double(tmp{1,1}{12,1})];
        end
        ggatime=tmp{1,1}{2,1};
        cursecs=round((datenum(ggatime,'HHMMSS')...
        - floor(datenum(ggatime,'HHMMSS')))*86400)+1;
        if exist('cursecsold','var')~=0 && exist('dt','var')==0
        dt=cursecs-cursecsold;
        disp(['frequency of data is ',num2str(dt),' s'])
        end
        cursecsold=cursecs;
    elseif strcmp(tline(4:6),'GSV')==1
        if numel(tline)-5>8
        tlinet=tline(8:end-5);
        values=textscan(tlinet,'%f','delimiter',',');
        values=values{:};
        sats=values(4:4:end-3);
        if strcmp(tline(3),'L')
            sats=sats-32; % + 32-64
        elseif strcmp(tline(3),'A')
            sats=sats+32+24;
        end
        if numel(sats)>0
        for ss=1:size(sats,1)
            if ~isnan(sats(ss)) && ~isnan(values(ss*4+3)) && exist('cursecs')~=0
                s1data(cursecs,sats(ss),1)=sats(ss);
                s1data(cursecs,sats(ss),2)=values(ss*4+1); % elv
                s1data(cursecs,sats(ss),3)=values(ss*4+2); % azi
                s1data(cursecs,sats(ss),4)=cursecs;
                s1data(cursecs,sats(ss),5)=values(ss*4+3); % S1
            end
        end
        end
        end
    end
end
if feof(fid)
    donefile=1;
end

if cnt==1
    disp('empty file')
    disp('try next one')
    return
end

if numel(sp3dir)>0
lattt=nanmedian(lat);
lontt=nanmedian(lon);
heighttt=nanmedian(height);
[staxyz(1),staxyz(2),staxyz(3)]=lla2ecef(lattt*pi/180,lontt*pi/180,heighttt);
clear lat lon height

curdt=datetime(datevecn,'convertfrom','datenum');
strday=char(datetime(curdt,'format','DDD'));
stryrl=char(datetime(curdt,'format','yyyy'));
[gpsw,sow]=datenum2gps(datevecn);
dow=sow/86400;

if sp3option==1
% OPTION 1
sp3str=[sp3dir,'/com/com',num2str(gpsw),num2str(round(dow)),'.sp3'];
elseif sp3option==2
% OPTION 2
sp3str=[sp3dir,'/igs/igs',num2str(gpsw),num2str(round(dow)),'.sp3'];
elseif sp3option==3
% OPTION 3
sp3str=[sp3dir,'/COD/COD0MGXFIN_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
elseif sp3option==4
% OPTION 4
sp3str=[sp3dir,'/GFZ/GFZ0MGXRAP_',stryrl,strday,'0000_01D_05M_ORB.SP3'];
end
 
disp('getting orbit info')
 
[txyz,xyz] = readsp3file(sp3str);
txyzsecs=txyz-txyz(1);
txyzsecs=txyzsecs.*86400;

snr_data=[];
for satind=1:92
    
s1datat=squeeze(s1data(:,satind,5));
secs=find(~isnan(s1datat(:)));
secs=secs-1; % because 1 is 0
if tempres~=0
    inspl=mod(secs,tempres)==0;
    secs=secs(inspl);
end
s1datat=s1datat(secs+1);

clear xyzt
if sum(~isnan(squeeze(xyz(satind,:,1))))>1 && numel(s1datat)>0
xyzt(:,1)=spline(txyzsecs,squeeze(xyz(satind,:,1)),secs);
xyzt(:,2)=spline(txyzsecs,squeeze(xyz(satind,:,2)),secs);
xyzt(:,3)=spline(txyzsecs,squeeze(xyz(satind,:,3)),secs);
else
    continue
end

ind=0;
snrdatat=[];
for tt=1:numel(secs)
[azi,elv]=gnss2azelv(staxyz,xyzt(tt,:),lattt,lontt);
if (azi>azilims(1) && azi<azilims(2)) && (elv>elvlims(1) && elv<elvlims(2))
ind=ind+1;
snrdatat(ind,1)=satind;
snrdatat(ind,2)=elv;
snrdatat(ind,3)=azi;
snrdatat(ind,4)=secs(tt);
snrdatat(ind,5)=NaN;
snrdatat(ind,6)=NaN;
snrdatat(ind,7)=s1datat(tt);
snrdatat(ind,8)=NaN;
snrdatat(ind,9)=datevecn+secs(tt)/86400;
end
end

snr_data=[snr_data;snrdatat];

end

clear snrdatat

else
    snr_data=[];
for ij=1:92
    inds=~isnan(s1data(:,ij,5));
    if sum(inds)>1
    s1datat=squeeze(s1data(inds,ij,:));
    inds=s1datat(:,2)>elvlims(1) & s1datat(:,2)<elvlims(2);
    s1datat=s1datat(inds,:);
    inds=s1datat(:,3)>azilims(1) & s1datat(:,3)<azilims(2);
    s1datat=s1datat(inds,:);
    % now fit 2nd order poly
    p2=polyfit(s1datat(:,4),s1datat(:,2),4);
    s1datat(:,2)=polyval(p2,s1datat(:,4));
    %newelv=polyval(p2,s1datat(:,4));
    %plot(s1datat(:,4),s1datat(:,2))
    %hold on
    %plot(s1datat(:,4),newelv)
    %return
    snr_data=[snr_data;s1datat];
    end
end
snr_data(:,7)=snr_data(:,5);
snr_data(:,5)=NaN;
snr_data(:,6)=NaN;
snr_data(:,8)=NaN;
end

if size(snr_data,1)>0
snr_data=sortrows(snr_data,4);
if exist(snrdir)==0
    mkdir(snrdir)
end
if exist([snrdir,'/',num2str(datevecn),'.mat'])~=0
    snr_datanew=snr_data;
    load([snrdir,'/',num2str(datevecn),'.mat'])
    snr_data=[snr_data;snr_datanew];
    snr_data(:,8)=[];
    snr_data(:,6)=[];
    snr_data(:,5)=[];
    snr_data=unique(snr_data,'rows');
    snr_data(:,7)=snr_data(:,5);
    snr_data(:,8)=NaN;
    snr_data(:,6)=NaN;
    snr_data(:,5)=NaN;
    snr_data=sortrows(snr_data,4);
end
lat=lattt;
lon=lontt;
height=heighttt;
save([snrdir,'/',num2str(datevecn),'.mat'],'snr_data','staxyz','lat','lon','height')
clear lat lon height
else
    disp('no data')
    return
end
clear snr_data snrdatanew

disp('** saved one day of data **')

end

end


