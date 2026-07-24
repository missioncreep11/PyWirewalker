function combine_cutoff(WWmeta,num,k)
% function to combine the profiles separated into consecutive files
% Bofu Zheng

dd0 = dir([WWmeta.aqdpath '*.mat']);
inds = 1:num:length(dd0);  % starting index for each group
inde = inds(2:end)-1;
inde = [inde length(dd0)]; % ending index for each group

%%
% process upcast
for q = 1:length(inds)-1

filename_1 = ['Profiles_upcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file1
filename_2 = ['Profiles_upcast_',WWmeta.name_aqd,'_',num2str(inds(q+1)),'_',num2str(inde(q+1)),'.mat']; % Name of the processed file2

aqd1 = load([WWmeta.propath_rearrange filename_1]);
aqd2 = load([WWmeta.propath_rearrange filename_2]);

pres_1 = aqd1.AQDprofiles_up{end}.Burst_Pressure;
pres_2 = aqd2.AQDprofiles_up{1}.Burst_Pressure;
if isfield(aqd1.AQDprofiles_up{end},'Burst_Time')
    time_1 = aqd1.AQDprofiles_up{end}.Burst_Time;
elseif isfield(aqd1.AQDprofiles_up{end},'Burst_MatlabTimeStamp')
    time_1 = aqd1.AQDprofiles_up{end}.Burst_MatlabTimeStamp;
else
    disp("Error: Could not find a time variable")
end
if isfield(aqd1.AQDprofiles_up{1},'Burst_Time')
    time_2 = aqd1.AQDprofiles_up{1}.Burst_Time;
elseif isfield(aqd1.AQDprofiles_up{1},'Burst_MatlabTimeStamp')
    time_2 = aqd1.AQDprofiles_up{1}.Burst_MatlabTimeStamp;
else
    disp("Error: Could not find a time variable")
end

if abs(pres_1(end)-pres_2(1))<10 & abs(time_1(end)-time_2(1))*86400<20% this means profile can be combined together
    disp(['start combination upcast; pres1:',num2str(pres_1(end)),', pres2:',num2str(pres_2(1))...
        ', file location: ',num2str(inds(q)),'-',num2str(inde(q))])
    tdata1=aqd1.AQDprofiles_up{end}.Burst_Time;
    tdata2=aqd2.AQDprofiles_up{1}.Burst_Time;
    fields=fieldnames(aqd1.AQDprofiles_up{end}); % identify fields
    for f=1:length(fields)   % loop with fields
        wh_field=fields{f};    % field name
        if (length(tdata1)==length(aqd1.AQDprofiles_up{end}.(wh_field)))
            if ~isfield(aqd1.AQDprofiles_up{end},wh_field)
                aqd1.AQDprofiles_up{end}.(wh_field)=[nan(size(aqd1.AQDprofiles_up{end}.Burst_Time,1),size(aqd2.AQDprofiles_up{1}.(wh_field),2));aqd2.AQDprofiles_up{1}.(wh_field)];
            elseif ~isfield(aqd2.AQDprofiles_up{end},wh_field)
                aqd1.AQDprofiles_up{end}.(wh_field)=[aqd1.AQDprofiles_up{end}.(wh_field);nan(size(aqd1.AQDprofiles_up{end}.Burst_Time,1),size(aqd1.AQDprofiles_up{end}.(wh_field),2))];
            elseif ~isfield(aqd2.AQDprofiles_up{end},wh_field) & ~isfield(aqd2.AQDprofiles_up{end},wh_field)
                aqd1.AQDprofiles_up{end}.(wh_field) = nan;
            else
                aqd1.AQDprofiles_up{end}.(wh_field)=[aqd1.AQDprofiles_up{end}.(wh_field);aqd2.AQDprofiles_up{1}.(wh_field)];
            end
        end
    end

    aqd2.AQDprofiles_up = aqd2.AQDprofiles_up(2:end);
    AQDprofiles_up = aqd1.AQDprofiles_up;
    save([WWmeta.propath_rearrange filename_1],'AQDprofiles_up', '-v7.3')
    AQDprofiles_up = aqd2.AQDprofiles_up;
    save([WWmeta.propath_rearrange filename_2],'AQDprofiles_up', '-v7.3')

end
end

%%
if k ==1  % process this when k==1
% process downcast
for q = 1:length(inds)-1
filename_1 = ['Profiles_downcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file1
filename_2 = ['Profiles_downcast_',WWmeta.name_aqd,'_',num2str(inds(q+1)),'_',num2str(inde(q+1)),'.mat']; % Name of the processed file2
% propath='/Users/warrbob/Desktop/Profile/'; % path to the save processed data - profile
aqd1 = load([WWmeta.propath_rearrange filename_1]);
aqd2 = load([WWmeta.propath_rearrange filename_2]);

pres_1 = aqd1.AQDprofiles_down{end}.Burst_Pressure;
pres_2 = aqd2.AQDprofiles_down{1}.Burst_Pressure;
if isfield(aqd1.AQDprofiles_up{end},'Burst_Time')
    time_1 = aqd1.AQDprofiles_up{end}.Burst_Time;
elseif isfield(aqd1.AQDprofiles_up{end},'Burst_MatlabTimeStamp')
    time_1 = aqd1.AQDprofiles_up{end}.Burst_MatlabTimeStamp;
else
    disp("Error: Could not find a time variable")
end
if isfield(aqd1.AQDprofiles_up{1},'Burst_Time')
    time_2 = aqd1.AQDprofiles_up{1}.Burst_Time;
elseif isfield(aqd1.AQDprofiles_up{1},'Burst_MatlabTimeStamp')
    time_2 = aqd1.AQDprofiles_up{1}.Burst_MatlabTimeStamp;
else
    disp("Error: Could not find a time variable")
end


if abs(pres_1(end)-pres_2(1))<10 & abs(time_1(end)-time_2(1))*86400<20% this means profile can be combined together
    disp(['start combination downcast; pres1:',num2str(pres_1(end)),', pres2:',num2str(pres_2(1))...
        ', file location: ',num2str(inds(q)),'-',num2str(inde(q))])
    tdata1=aqd1.AQDprofiles_down{end}.Burst_Time;
    tdata2=aqd2.AQDprofiles_down{1}.Burst_Time;
    fields=fieldnames(aqd1.AQDprofiles_down{end}); % identify fields
    for f=1:length(fields)   % loop with fields
        wh_field=fields{f};    % field name
        if (length(tdata1)==length(aqd1.AQDprofiles_down{end}.(wh_field)))
            aqd1.AQDprofiles_down{end}.(wh_field)=[aqd1.AQDprofiles_down{end}.(wh_field);aqd2.AQDprofiles_down{1}.(wh_field)];
        end
    end

    aqd2.AQDprofiles_down = aqd2.AQDprofiles_down(2:end);
    AQDprofiles_down = aqd1.AQDprofiles_down;
    save([WWmeta.propath_rearrange filename_1],'AQDprofiles_down', '-v7.3')
    AQDprofiles_down = aqd2.AQDprofiles_down;
    save([WWmeta.propath_rearrange filename_2],'AQDprofiles_down', '-v7.3')
end
end
end
    
end