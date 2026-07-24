% This function cuts out sections of a continous ADCP record to simulate
% duty cyling the instrument. Currently set up to simulate the 1/3rd duty
% cycle used on SWOT moorings

function merge_signature_SWOT_emulator(WWmeta,q,num)
% q - starting position
% num - number of file for combination
% 
% Devon Northcott, Sep. 2023
% Based on previous version by Arnaud LeBoyer and Bofu

beg=zeros(1,num);
cell_Data=struct([]);
for l=1:num
    filename = WWmeta.sortedname{q+l-1};
    load([WWmeta.aqdpath filename]);
    beg(l)=Data.Burst_Time(1);
    cell_Data{l}=Data;
    cell_Config{l}=Config;
end
[~,I]=sort(beg);
Fields=fields(Data);
AllData=[cell_Data{I}];

% input to emulate SWOT 1/3rd duty cycle
Duty_cycle_mask = zeros(size(Data.Burst_Time));
splits = Data.Burst_Time(1):30/1440:Data.Burst_Time(end); % Set to adjust "on time"
ind = ceil(rand(1)*3);
for n = 1:length(splits)-1
    mask = Data.Burst_Time>splits(n) & Data.Burst_Time<splits(n+1);
    ind=ind+1;
    if mod(ind,3)==0
        Duty_cycle_mask(mask)=1;
    end
end
Duty_cycle_mask = Duty_cycle_mask==1;

AllData1=struct();
for f=1:length(Fields)
    field=Fields{f};
    AllData1.(field)=vertcat(AllData(:).(field)(Duty_cycle_mask,:));
end

name_file = [WWmeta.name_aqd,'_',num2str(q),'_',num2str(q+num-1)];
eval([name_file '=AllData1;']);
save([WWmeta.matpath name_file '.mat'],name_file, '-v7.3')
end