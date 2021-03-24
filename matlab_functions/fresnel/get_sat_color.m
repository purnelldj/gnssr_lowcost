function col = get_sat_color(sat)
%function col = get_sat_color(sat)
% GNSSH2O
% Author: Kristine Larson
% This function loads the RGB color corresponding to satellite (sat) and
% normalizes it into [0-1,0-1,0-1] values for matlab
% kristine wrote this
% default is black 
col = [0 0 0];
colormatrix=load([getenv('INSTRUCTIONS') 'satellite_colors.txt']);
i=find(sat == colormatrix(:,1));
if length(i) > 0
  col=colormatrix(i,2:4)/255; % normalize RGB values to 0-1
end
  
end