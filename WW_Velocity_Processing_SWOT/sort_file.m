function WWmeta = sort_file(WWmeta)
% to sort the raw .mat file in the right order in time
%
% Bofu Zheng, Nov. 23, 2020

dd0 = dir([WWmeta.aqdpath '*.mat']);  % info of all .mat files
% load time and sort files
time = [];
name = [];  % name of all files, help with future loading
for i = 1:length(dd0)
    filename = dd0(i).name;
    load([WWmeta.aqdpath,filename]);
    if isfield(Data,'Burst_Time')
        time = [time;Data.Burst_Time(1)];
    elseif isfield(Data,'Burst_MatlabTimeStamp')
        time = [time;Data.Burst_MatlabTimeStamp(1)];
    else
        disp("Error: Could not find a time variable")
    end
    name{i} = dd0(i).name; 
end

% sort time
[~,I]=sort(time);
names = name(I);  % sorted name
WWmeta.sortedname = names;
end