clear all; close all;

statsdir='~/data/rv3sa/snr';

files=dir(statsdir);

heights=[];
lats=[];
lons=[];
staxyzs=[];
names={files.name};
for ii=1:size(names,2)
    tmpfile=char(names(ii));
    if numel(tmpfile)>4
    if strcmp(tmpfile(end-3:end),'.mat')
        load([statsdir,'/',tmpfile])
        lats=[lats lat];
        lons=[lons lon];
        heights=[heights height];
        staxyzs=[staxyzs; staxyz];
    end
    end
end

format long
disp(['lat is ',num2str(nanmedian(lats))]);
disp(['lon is ',num2str(nanmedian(lons))]);
disp(['height is ',num2str(nanmedian(heights))]);
disp(['staxyz is ',num2str(nanmedian(staxyz,1))]);
