function [azimuth,elev] = gnss2azelv(staXYZ,SatPos,lat,lon)

% converts station position and satellite position to azimuth and elevation
% angles
% written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost

dx=SatPos(1)-staXYZ(1);
dy=SatPos(2)-staXYZ(2);
dz=SatPos(3)-staXYZ(3);
%lla=ecef2lla([ecestaXYZ(1) staXYZ(2) staXYZ(3)]);
%range=sqrt(dx^2+dy^2+dz^2);

ll=[lat;lon];

trans_matrix = zeros(3, 3);
trans_matrix(1,1)=-sind(ll(2));
trans_matrix(2,1)=-sind(ll(1))*cosd(ll(2));
trans_matrix(3,1)=cosd(ll(1))*cosd(ll(2));
trans_matrix(1,2)=cosd(ll(2));
trans_matrix(2,2)=-sind(ll(1))*sind(ll(2));
trans_matrix(3,2)=cosd(ll(1))*sind(ll(2));
trans_matrix(1,3)=0;
trans_matrix(2,3)= cosd(ll(1));
trans_matrix(3,3)= sind(ll(1));

obsvec=[dx;dy;dz];

rss = trans_matrix*obsvec;

azimuth = atan2(rss(1, :), rss(2, :))*180/pi;
if azimuth<0
    azimuth=360+azimuth;
end
elev = asind(rss(3, :)/sqrt(sum((rss).^2)));

end





