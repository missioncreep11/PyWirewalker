function [up,down,dataup,datadown] = get_aqd_2G(data,thhold)
% modified by bzheng
% thhold: threshold value to determine whether a section is too short to be
% a profile
% 
% rename variable to make it easy
pdata=double(data.Burst_Pressure);
if isfield(data,'Burst_Time')
    tdata=data.Burst_Time;
elseif isfield(data,'Burst_MatlabTimeStamp')
    tdata=data.Burst_MatlabTimeStamp;
else
    disp("Error: Could not find a time variable")
end
dtdata = [0;diff(tdata)];

% buid a filter 
dt=median(diff(tdata)); % sampling period
T=tdata(end)-tdata(1);  % length of the record
disp('check if time series is shorter than 3 hours')
if T<3/24  
    warning('time serie is less than 3 hours, very short for data processing, watch out the results')
end

disp('smooth the pressure to define up and down cast')
Nb  = 3; % filter order
fnb = 1/(2*dt); % Nyquist frequency
fc  = 1/200/dt; % 600 dt (give "large scale patern") 
[b,a]= butter(Nb,fc/fnb,'low');
filt_pdata=filtfilt(b,a,pdata);
dfilt_pdata=diff(filt_pdata);
pddfilt_pdata = dfilt_pdata(1:end-1).*dfilt_pdata(2:end);
pind = find(pddfilt_pdata <= 0);
csind = [1; pind(1:end)+1];     % long column of cast start time
ceind = [pind(1:end); length(pdata)];   % cast endings
thind = [];
for i = 1:length(csind)
    if ceind(i)-csind(i) < thhold  % not a profile
        thind = [thind i]; 
    end
end
csind(thind) = [];
ceind(thind) = [];
fprintf('identify %i profile \n',round(length(csind)/2))
%%
% Added by Devon Northcott, 9/23, to account for big jumps in time,
% especally those caused by duty cycling the ADCP. Starts a new profile if
% more than 30s (~15m) of data is missing.
time_jump_mask = find(dtdata>30/86400);
csind = [csind;time_jump_mask+1];
ceind = [ceind;time_jump_mask-1];
csind = unique(csind);ceind = unique(ceind);

slope=nan*pdata;
for i=1:length(csind)
    [~,mI]=min(pdata(csind(i):ceind(i)));
    [~,MI]=max(pdata(csind(i):ceind(i)));
    if MI>mI
        slope(csind(i):ceind(i))=0;%downcast
    else
        slope(csind(i):ceind(i))=1;%upcast
    end
end
down = find(slope==0);
up = find(slope==1);


% Plot it
figure(1);clf;
plot(tdata(down),pdata(down),'b.')
hold on
plot(tdata(up),pdata(up),'r.')
set(gca,'ydir','reverse')
title('Verify rising/falling data separation');
ylabel('depth (m)')
set(gca,'xtick',tdata(1):(tdata(end)-tdata(1))/4:tdata(end),'tickdir','in');
datetick('x','dd HH:MM','keepticks');
xlabel('time (DD HH:MM)')
xlim([tdata(1) tdata(end)])

% Added by Devon Northcott, 7/24, Cluge to fix situations where HR mode
% data has one more or fewer data points than non-HR data.
fields=fieldnames(data);
if (sum(strcmp(fields,'IBurstHR_Time'))==1 | sum(strcmp(fields,'IBurstHR_MatlabTimeStamp'))==1)
    if isfield(data,'IBurstHR_Time')
        HRtdata=data.IBurstHR_Time;
    elseif isfield(data,'IBurstHR_MatlabTimeStamp')
        HRtdata=data.IBurstHR_MatlabTimeStamp;
    else
        disp("Error: Could not find a time variable")
    end
 
    if length(tdata)~=length(HRtdata)
        HRfields_cell = strfind(fields,'IBurstHR');
        HRfields_bool = zeros(size(HRfields_cell));
        for i = 1:length(HRfields_cell)
            HRfields_bool(i) = ~isempty(HRfields_cell{i});
        end
        HRfields = fields;
        HRfields(~HRfields_bool)=[];
        len_diff = length(tdata)-length(HRtdata);
        if len_diff<0
            for f=1:length(HRfields)
                wh_field=HRfields{f};
                data.(wh_field) = data.(wh_field)(1:end+len_diff,:);
            end
        elseif len_diff>0
            for f=1:length(HRfields)
                wh_field=HRfields{f};
                data.(wh_field) = [data.(wh_field);nan(len_diff,size(data.(wh_field),2))];
            end
        end
    end
end

% if nargout==4
    for f=1:length(fields)
        wh_field=fields{f};
        if (length(tdata)==length(data.(wh_field))||(length(tdata)+1)==length(data.(wh_field)))
            switch length(size(data.(wh_field)))
                case 2
                    dataup.(wh_field)=data.(wh_field)(up,:);
                    datadown.(wh_field)=data.(wh_field)(down,:);
                case 3
                    dataup.(wh_field)=data.(wh_field)(up,:,:);
                    datadown.(wh_field)=data.(wh_field)(down,:,:);
                case 4
                    dataup.(wh_field)=data.(wh_field)(up,:,:,:);
                    datadown.(wh_field)=data.(wh_field)(down,:,:,:);
            end
        end
    end
% end

end

