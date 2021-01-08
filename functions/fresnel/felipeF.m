function [firstF] = felipeF(freq, e, h, theta)
% input elevation angle (degrees)
%       h (antenna height in meters)
%       theta (azimuth in degrees)
% satellite elevation angle (assumes horizontal, untilted surface)
%

%SOME GPSCONSTANTS	
CLIGHT = 299792458;             % speed of light, m/sec
FREQ = [1575.42e6; 1227.6e6];   %GPS frequencies, Hz
CYCLE = CLIGHT./FREQ;           %wavelength per cycle (m/cycle)
RAD2M = CYCLE/2/pi;             %(m)

% ------------------
% ------------------
% delta = locus of points corresponding to a fixed delay;
% typically the first Fresnel zone is is the 
% "zone for which the differential phase change across
% the surface is constrained to lambda/2" (i.e. 1/2 the wavelength)
delta = CYCLE(freq)/2; 	% [meters]


% equations from Katzberg 1996 paper
%q = sqrt(2.*delta.*h.*sind(e));	% [meters]
%A = q./(sind(e).^2);
%B = q./sind(e);

% Penny Axelrad's approximations (Equation 15 in Raytheon report)
%A = sqrt(CYCLE(freq)*h./(sind(e).^3)) 
%B = sqrt(CYCLE(freq)*h./(sind(e))) 


% Felipe Nievinski values, from Larson and Nievinski, 2013
sin_elev = sind(e);
%d = n .* wavelength ./ 2;
d = delta; % 
B = sqrt( (2 .* d .* h ./ sin_elev) + (d ./ sin_elev).^2 ) ;
A = B ./ sin_elev ;


% determine distance to ellipse center and plot
% Katzberg, Equation 6:
center(:,1) = (h + delta./sind(e))./tand(e);
% Penny Axelrad, Equation 12:
center(:,2) = (delta + h*sind(e))./(sind(e).*tand(e));
%center(:,3) = h./tand(e);

[firstF]=[A',B',center(:,1)];

end