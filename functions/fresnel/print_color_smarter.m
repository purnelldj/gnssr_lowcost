function print_color_smarter(fid, PRN)
% function print_color(fid, PRN)
% called by google earth code
% it would be nice if someone made this less stupid. KL
% blue green red, whereas matlab uses red green blue
% no one bothered to provide a link for htis - maybe this will help
% http://msdn.microsoft.com/en-us/library/system.drawing.color.aspx
% should be else if but life is short
col3=get_sat_color(PRN);
if sum(col3) == 0
  %disp('print nothing')
else
  newval = ge_color(col3);
  fprintf(fid,['FF' newval]);
end
end