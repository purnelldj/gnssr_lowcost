function [txyz,xyz] = readsp3file(sp3str)

% this function reads an sp3 orbit file and outputs timex and x,y,z
% positions of satellites
% written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost

fid=fopen(sp3str);

tline=fgets(fid);
while ~strcmp(tline(1),'*')
    tline=fgets(fid);
end

txyz=NaN(86400/(5*60)+1,1); % 5 mins would be max so preallocating for speed
xyz=NaN(92,numel(txyz),3);
tt=0;
while ~feof(fid)
    
    tt=tt+1;
    txyz(tt)=datenum(str2double(tline(4:7)),str2double(tline(9:10)),str2double(tline(12:13)),...
        str2double(tline(15:16)),str2double(tline(18:19)),str2double(tline(21:31)));
    tline=fgets(fid);
    while strcmp(tline(1),'*')==0 && ~feof(fid)
        if strcmp(tline(2),'G')
            satidt=str2double(tline(3:4));
        elseif strcmp(tline(2),'R')
            satidt=str2double(tline(3:4))+32;
        elseif strcmp(tline(2),'E')
            satidt=str2double(tline(3:4))+56;
        else
            tline=fgets(fid);
            continue
        end
        xyz(satidt,tt,:)=[str2double(tline(5:18)) str2double(tline(19:32)) str2double(tline(33:46))];
        tline=fgets(fid);
    end
end

if tt~=numel(txyz)
    txyz(tt+1:end)=[];
    xyz(:,tt+1:end,:)=[];
end

if mod(txyz(1),1)~=0
    disp('starting point not start of day - bad')
end

xyz=xyz.*1000; % convert from km to m - IMPORTANT

fclose(fid);

end

