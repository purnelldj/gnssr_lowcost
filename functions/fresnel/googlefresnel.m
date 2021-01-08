function googlefresnel(lat,lon,el,az,PRN,savedir, ht, freq)
%function googlefresnel(lat,lon,el,az,PRN,station, ht, freq)
% inputs
% latitude (deg)
% longitude (deg)
% el is elevation angle (deg)
% az is azimuth angle (deg)
% PRN is satellite number
% savedir is directory to save to
%Radius of Earth, average
R=6378.14; %km
% define quadrants
quadlist = ['ne'; 'nw'; 'sw'; 'se'];
q = 1;
if az < 90, q= 1; end;
if az >= 90 & az < 180, q= 4; end;
if az >= 180 & az < 270 , q= 3; end;
if az > 270, q = 2; end;
% you will need to change the output location for the files
dirdir = [savedir,'/'];
%if (~exist(dirdir))
%  unix(['mkdir ' dirdir]);
%  unix(['chmod g+rw ' dirdir]);
%end

F = felipeF(freq, el, ht, az);
a=F(1);
b=F(2);
center=F(3);
%Calculate the relative x and y positions for points along the Fresnel
%ellipse.  
%Convert az to proper angle for use in ellipse. azimuth is typically
%measured clockwise with north as zero, ellipse is a cartesian system which
%is counter-clockwise and east is zero.

azcart=360-az+90;
if azcart>360
    azcart=azcart-360;
end
%convert aznew to radians
azcart=azcart*pi/180;

[x y]=ellipseGE(a,b,azcart,center*cos(azcart),center*sin(azcart));

%for each of the x-y coordinates, calculate the distance from the antenna,
%the bearing angle relative to the antenna in order to solve for the
%latitude and longitude of each point.


d=sqrt(x.^2+y.^2); %meters ; 
d=d./1000; %km

%Calculate bearing angle. This is reference similiarly as azimuth i.e.
%clockwise from north.

theta=atan2(x,y);
k=find(theta<0);
theta(k)=theta(k)+2*pi;
theta=theta*180/pi;

%new lat and lon
latnew=asin(sind(lat).*cos(d./R)+cosd(lat).*sin(d./R).*cosd(theta));
lonnew=lon+180./pi.*(atan2(sind(theta).*sin(d./R).*cosd(lat),cos(d./R)-sind(lat).*sin(latnew)));
latnew=latnew.*180./pi;
data=[lonnew;latnew];

azd=azcart*180/pi;
filename= ['h' num2str(ht) '_' quadlist(q,1:2) num2str(PRN) '_e', num2str(el) '_a' num2str(round(az)),'.kml'] ;
writefresnel_ge(filename, dirdir, PRN, az, el,data)

end