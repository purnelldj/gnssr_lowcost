function [gpsweek,sow] = datenum2gps(datenm)

% written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost

curdt = datetime(datenm,'convertfrom','datenum');
jd = juliandate(curdt);
curdtgps = datetime(1980,1,6);
jdgps = juliandate(curdtgps);
nweek = fix((jd-jdgps)/7);
sow = (jd - (jdgps+nweek*7)) * 3600*24;
% rollover = fix(nweek/1024);  % rollover every 1024 weeks
gpsweek = nweek;


end