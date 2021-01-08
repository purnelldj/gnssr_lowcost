function writefresnel_ge(filename, dirdir, PRN, az, el,data)
%function writefresnel_ge(filename, dirdir, PRN, az, el,data)
% inputs filename, output directory, PRN number, azimuth, elevation angles
% ellipse coordinates are in the variable data
%disp(['opening ' dirdir filename])
fid = fopen([dirdir filename], 'w');
fprintf(fid,'<Folder>\n<name>');
fprintf(fid,'SV %2.0f %3f %2f', PRN, az, el); %SV ## AZd Eld
fprintf(fid,'</name>\n<visibility>1</visibility>\n<Placemark>\n<name>circle</name>\n<visibility>1</visibility>\n<Style>\n<geomColor>');
print_color_smarter(fid, PRN);

fprintf(fid, '</geomColor>\n<geomScale>2</geomScale></Style>\n<LineString>\n<coordinates>\n');
% write out the file - stored very strangely
[nr,nc]=size(data);
for i=1:nc
  fprintf(fid, '%10.7f,%10.7f,0 \n', data(1,i), data(2,i));
end
fprintf(fid,'</coordinates>\n</LineString>\n</Placemark>\n</Folder>');
fclose(fid);

end