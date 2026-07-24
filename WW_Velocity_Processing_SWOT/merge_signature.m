function merge_signature(WWmeta,q,num)
% q - starting position
% num - number of file for combination
% 
% Bofu Zheng, Nov. 2020
% Based on previous version by Arnaud LeBoyer

beg=zeros(1,num);
cell_Data=struct([]);
for l=1:num
    filename = WWmeta.sortedname{q+l-1};
    load([WWmeta.aqdpath filename]);
    if isfield(Data,'Burst_Time')
        beg(l) = Data.Burst_Time(1);
    elseif isfield(Data,'Burst_MatlabTimeStamp')
        beg(l) = Data.Burst_MatlabTimeStamp(1);
    else
        disp("Error: Could not find a time variable")
    end
    cell_Data{l}=Data;
    cell_Config{l}=Config;
end
[~,I]=sort(beg);
Fields=fields(Data);
AllData=[cell_Data{I}];

AllData1=struct();
for f=1:length(Fields)
    field=Fields{f};
    AllData1.(field)=vertcat(AllData(:).(field));
end

name_file = [WWmeta.name_aqd,'_',num2str(q),'_',num2str(q+num-1)];
eval([name_file '=AllData1;']);
save([WWmeta.matpath name_file '.mat'],name_file, '-v7.3')
end