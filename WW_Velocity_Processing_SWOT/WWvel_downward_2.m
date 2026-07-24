function out = WWvel_downward_2(WWmeta,variables,k)
% out{1} = [];  % E-W velocity
% out{2} = [];  % N-s velocity
% out{3} = [];  % time
% out{4} = [];  % depth
% bofu zheng

%% set variables
if nargin<2
    num      = 20;   % 20 files in a group 
    blockdis = 0.1;  % blocking distance 0.1 m
    cellsize = 0.25; % cell size 0.25 m
    saprate  = 16;   % sampling rate 16 Hz
    boxs  = 0.25; % box size 0.25 m
    z_max    = 100;  % profile depth 100 m
    direction = 'down';
else
    num      = variables.NUM_combining_files;
    blockdis = variables.blockdis;
    cellsize = variables.cellsize;
    saprate  = variables.saprate;
    boxs  = variables.boxsize;
    z_max    = variables.z_max;
    direction = variables.direction;
end

%%

out{1} = [];  % E-W velocity
out{2} = [];  % N-s velocity
out{3} = [];  % time
out{4} = [];  % depth
out{5} = [];  % E-W shear
out{6} = [];  % N-S shear

dd0 = dir([WWmeta.aqdpath '*.mat']);
inds = 1:num:length(dd0);  % starting index for each group
inde = inds(2:end)-1;
inde = [inde length(dd0)]; % ending index for each group

%% analysis for upcast

vele_up = [];  % final product, E/W velocity
veln_up = [];  % final product, N/S velocity
time_up = [];  % final product, time
shearu_up = [];
shearv_up = [];
surf_vel = [];

for q = 1:length(inds)
    filename = ['Profiles_upcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file
    load([WWmeta.propath_rearrange filename]);
    % find long enough profile
%     AQDprofiles_up = AQDprofiles;
    index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_up))>0.4*z_max);
    
    %% load data
    for k = 1:length(index)
        time = AQDprofiles_up{index(k)}.Burst_Time;
        pres = (AQDprofiles_up{index(k)}.Burst_Pressure);%/10*100000)/(9.8*1025);
        temp = (AQDprofiles_up{index(k)}.Burst_Temperature);
        rho = gsw_rho(ones(size(temp))*33,temp,pres);
        dpth = (AQDprofiles_up{index(k)}.Burst_Pressure/10*100000)./(9.8*rho);
        
        if isfield(AQDprofiles_up{index(k)},'Burst_CorBeam1')
            corr1 = double(AQDprofiles_up{index(k)}.Burst_CorBeam1);
            corr2 = double(AQDprofiles_up{index(k)}.Burst_CorBeam2);
            corr3 = double(AQDprofiles_up{index(k)}.Burst_CorBeam3);
            corr4 = double(AQDprofiles_up{index(k)}.Burst_CorBeam4);
        else
            corr1 = double(squeeze(AQDprofiles_up{index(k)}.Burst_Correlation_Beam(:,1,:)));
            corr2 = double(squeeze(AQDprofiles_up{index(k)}.Burst_Correlation_Beam(:,2,:)));
            corr3 = double(squeeze(AQDprofiles_up{index(k)}.Burst_Correlation_Beam(:,3,:)));
            corr4 = double(squeeze(AQDprofiles_up{index(k)}.Burst_Correlation_Beam(:,4,:)));
        end
        
        beam_vel = [];
        if isfield(AQDprofiles_up{index(k)},'Burst_Velocity_Beam')
            beam_vel = AQDprofiles_up{index(k)}.Burst_Velocity_Beam;
        else
            beam_vel(:,1,:) = permute(AQDprofiles_up{index(k)}.Burst_VelBeam1,[1,3,2]);
            beam_vel(:,2,:) = permute(AQDprofiles_up{index(k)}.Burst_VelBeam2,[1,3,2]);
            beam_vel(:,3,:) = permute(AQDprofiles_up{index(k)}.Burst_VelBeam3,[1,3,2]);
            beam_vel(:,4,:) = permute(AQDprofiles_up{index(k)}.Burst_VelBeam4,[1,3,2]);
        end
        
        m = corr1<50 | corr2<50 | corr3<50 | corr4<50;

%         amp = double(squeeze(AQDprofiles_up{index(k)}.Burst_Amplitude_Beam(:,:,:)));
        
        pitch = AQDprofiles_up{index(k)}.Burst_Pitch;
        roll = AQDprofiles_up{index(k)}.Burst_Roll;
        heading = AQDprofiles_up{index(k)}.Burst_Heading;
        acce = AQDprofiles_up{index(k)}.Burst_Accelerometer;
        NCells = AQDprofiles_up{index(k)}.Burst_NCells(1); % number of cells
        Npings = length(time);   % number of pings
        
        phi = -[65,65,65,65]*pi/180;
        azi = [0,-90,180,90]*pi/180;

        clear bX bY bZ
        for ibeam=1:4
          [bX(:,ibeam),bY(:,ibeam),bZ(:,ibeam)] = GetUnitVectors(phi(ibeam),azi(ibeam),roll*pi/180,pitch*pi/180);  
        end

        ranges = (ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis)/cosd(25);
        z_coords = zeros(size(ranges,1),size(ranges,2),4);
        for n = 1:size(ranges,1)
            z_coords(n,:,:) = -dpth(n) - transpose(ranges(n,:)).*bZ(n,:);
        end
     
        % build dpth matrix
        if strcmp(direction,'down')
            dpth_temp = -dpth*ones(1,NCells)-ones(Npings,1)*cellsize*double(0:1:NCells-1) - blockdis;
        elseif strcmp(direction,'up')
            dpth_temp = -dpth*ones(1,NCells)+ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis;
        end
        dpth_temp(find(m == 1)) = nan;
        dpth_temp_all = dpth_temp(:);
        dpth_temp_shear = dpth_temp(:,2:end-1);
        dpth_temp_shear_all = dpth_temp_shear(:);
        
        beamshear = permute((diff(permute(beam_vel(:,:,1:end-1),[1,3,2]),1,2)-fliplr(diff(fliplr(permute(beam_vel(:,:,2:end),[1,3,2])),1,2)))...
            ./(diff(z_coords(:,1:end-1,:),1,2)-fliplr(diff(fliplr(z_coords(:,2:end,:)),1,2))),[1,3,2]);
        
        equal_depth_vel = zeros(size(beam_vel));
        for b = 1:4
            for n = 1:size(beam_vel,1)
                equal_depth_vel(n,b,:) = interp1(z_coords(n,:,b),squeeze(beam_vel(n,b,:)),dpth_temp(n,:));
            end
        end
        
        equal_depth_shear = zeros(size(beam_vel,1),size(beam_vel,2),size(beam_vel,3)-2);
        for b = 1:4
            for n = 1:size(beam_vel,1)
                equal_depth_shear(n,b,:) = interp1(z_coords(n,2:end-1,b),squeeze(beamshear(n,b,:)),dpth_temp_shear(n,:));
            end
        end
        
        equal_depth_ENU = Beam2ENU(squeeze(equal_depth_vel(:,1,:)),squeeze(equal_depth_vel(:,2,:)),...
            squeeze(equal_depth_vel(:,3,:)),squeeze(equal_depth_vel(:,4,:)), pitch, roll, heading);
        
        equal_depth_shear_ENU = Beam2ENU(squeeze(equal_depth_shear(:,1,:)),squeeze(equal_depth_shear(:,2,:)),...
            squeeze(equal_depth_shear(:,3,:)),squeeze(equal_depth_shear(:,4,:)), pitch, roll, heading);
        
        vele{1} = squeeze(equal_depth_ENU(:,1,:));
        vele{2} = squeeze(equal_depth_ENU(:,2,:));
        vele{3} = squeeze(equal_depth_ENU(:,3,:));

%         shearu1 = (diff(vele{1}(:,1:end-1),1,2)-fliplr(diff(fliplr(vele{1}(:,2:end)),1,2)))/(2*cellsize);
%         shearv1 = (diff(vele{2}(:,1:end-1),1,2)-fliplr(diff(fliplr(vele{2}(:,2:end)),1,2)))/(2*cellsize);
        shearu1 = smoothdata(smoothdata(squeeze(equal_depth_shear_ENU(:,1,:)),2,'movmean',1),'gaussian',16);
        shearv1 = smoothdata(smoothdata(squeeze(equal_depth_shear_ENU(:,2,:)),2,'movmean',1),'gaussian',16);
        shearu1(1:2,:)= NaN;
        shearv1(1:2,:)= NaN;
        shearu1(:,1:3)= NaN;
        shearv1(:,1:3)= NaN;
        
        % extract velocity that has low correlation
        vele{1}(find(m == 1)) = nan;
        vele{2}(find(m == 1)) = nan;
        vele{3}(find(m == 1)) = nan;
        
        % calculate ww motion velocity
        velemc = WWcorr(time,dpth,acce,pitch,roll,heading,vele,saprate);
        
        % build time matrix
        time_temp = time*ones(1,NCells);
        time_temp(find(m==1)) = nan;
        time_temp_all = time_temp(:);

        % all velocity 
        velemcE_temp_all = velemc{1}(:);  % E/W
        velemcN_temp_all = velemc{2}(:);  % N/S
        shearu1_temp = shearu1(:);
        shearv1_temp = shearv1(:);
        
        %% box averaging 
        z = (0:-boxs:-z_max)';
        temp_pro = nan(length(z),2);  % each profile, first column: E/W; second column: N/S
        temp_pro_shear = nan(length(z),2);
        temp_time   = nan(length(z),1);  % box-averaged time
        for i = 1:length(z)
            index_box = [];
            index_box = find(z(i)-boxs/2<dpth_temp_all & dpth_temp_all<=z(i)+boxs/2);
            if length(index_box)>0 & sum(~isnan(velemcE_temp_all(index_box)))>10 &...
                    sum(~isnan(velemcN_temp_all(index_box)))>10
                temp_pro(i,1)  = nanmean(velemcE_temp_all(index_box));
                temp_pro(i,2)  = nanmean(velemcN_temp_all(index_box));
                temp_time(i,1) = nanmean(time_temp_all(index_box));                       
            end
        end
        clear i
        
        for i = 1:length(z)
            index_box = [];
            index_box = find(z(i)-boxs/2<dpth_temp_shear_all & dpth_temp_shear_all<=z(i)+boxs/2);
            if length(index_box)>10 & sum(~isnan(shearu1_temp(index_box)))>10 &...
                    sum(~isnan(shearv1_temp(index_box)))>10
                temp_pro_shear(i,1)  = nanmean(shearu1_temp(index_box));
                temp_pro_shear(i,2)  = nanmean(shearv1_temp(index_box));                       
            end
        end
        clear i
        
        vele_up = [vele_up temp_pro(:,1)];
        veln_up = [veln_up temp_pro(:,2)];
        shearu_up = [shearu_up temp_pro_shear(:,1)];
        shearv_up = [shearv_up temp_pro_shear(:,2)];
        time_up = [time_up temp_time];
        
        disp([filename,' total #= ',num2str(size(time_up,2))]);
    end
end

%%
out{1} = vele_up;
out{2} = veln_up;
out{3} = time_up;
out{4} = -z*ones(1,size(time_up,2));
out{5} = shearu_up;
out{6} = shearv_up;
out{7} = surf_vel;

%% analysis for downcast - still needs work to be done
if k == 1

vele_down = [];
veln_down = [];
time_down = [];
for q = 1:length(inds)
    filename = ['Profiles_downcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file1
    load([WWmeta.propath_rearrange filename]);
    % find long enough profile
%     AQDprofiles_down = AQDprofiles;
    index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_down))>0.9*z_max);
    
    %% load data
    for k = 1:length(index)
        time = AQDprofiles_down{index(k)}.Burst_Time;
        dpth = AQDprofiles_down{index(k)}.Burst_Pressure;
        if isfield(AQDprofiles_down{index(k)},'Burst_CorBeam1')
            corr1 = double(AQDprofiles_down{index(k)}.Burst_CorBeam1);
            corr2 = double(AQDprofiles_down{index(k)}.Burst_CorBeam2);
            corr3 = double(AQDprofiles_down{index(k)}.Burst_CorBeam3);
            corr4 = double(AQDprofiles_down{index(k)}.Burst_CorBeam4);
        else
            corr1 = double(squeeze(AQDprofiles_down{index(k)}.Burst_Correlation_Beam(:,1,:)));
            corr2 = double(squeeze(AQDprofiles_down{index(k)}.Burst_Correlation_Beam(:,2,:)));
            corr3 = double(squeeze(AQDprofiles_down{index(k)}.Burst_Correlation_Beam(:,3,:)));
            corr4 = double(squeeze(AQDprofiles_down{index(k)}.Burst_Correlation_Beam(:,4,:)));
        end
        
        if isfield(AQDprofiles_down{index(k)},'Burst_VelEast')
            vele{1} = AQDprofiles_down{index(k)}.Burst_VelEast;
            vele{2} = AQDprofiles_down{index(k)}.Burst_VelNorth;
            vele{3} = 1/2*(AQDprofiles_down{index(k)}.Burst_VelUp1 + ...
                AQDprofiles_down{index(k)}.Burst_VelUp2);  % raw vel z
        else
            vele{1} = squeeze(AQDprofiles_down{index(k)}.Burst_Velocity_ENU(:,1,:));
            vele{2} = squeeze(AQDprofiles_down{index(k)}.Burst_Velocity_ENU(:,2,:));
            vele{3} = 1/2*(squeeze(AQDprofiles_down{index(k)}.Burst_Velocity_ENU(:,3,:)) + ...
                squeeze(AQDprofiles_down{index(k)}.Burst_Velocity_ENU(:,4,:)));
        end
        pitch = AQDprofiles_down{index(k)}.Burst_Pitch;
        roll = AQDprofiles_down{index(k)}.Burst_Roll;
        heading = AQDprofiles_down{index(k)}.Burst_Heading;
        acce = AQDprofiles_down{index(k)}.Burst_Accelerometer;
        NCells = AQDprofiles_down{index(k)}.Burst_NCells(1); % number of cells
        Npings = length(time);   % number of pings
    
        % extract velocity that has low correlation
        m = corr1<50 | corr2<50 | corr3<50 | corr4<50;
        vele{1}(find(m == 1)) = nan;
        vele{2}(find(m == 1)) = nan;
        vele{3}(find(m == 1)) = nan;
        
        vele{1}(:,1:5) = nan;
        vele{2}(:,1:5) = nan;
        vele{3}(:,1:5) = nan;
        
        % calculate ww motion velocity
        velemc = WWcorr(time,dpth,acce,pitch,roll,heading,vele,saprate);
        
%         ranged{1} = velemc{1} - nanmean(velemc{1},2)*ones(1,NCells);  % range dependent vel
%         ranged{2} = velemc{2} - nanmean(velemc{2},2)*ones(1,NCells);  % range dependent vel
%         ranged{3} = velemc{3} - nanmean(velemc{3},2)*ones(1,NCells);  % range dependent vel
%         
%         rangea{1} = nanmean(velemc{1},2)*ones(1,NCells);  % range averaged vel
%         rangea{2} = nanmean(velemc{2},2)*ones(1,NCells);  % range averaged vel
%         rangea{3} = nanmean(velemc{3},2)*ones(1,NCells);  % range averaged vel
            
        % build dpth matrix
        if strcmp(direction,'down')
            dpth_temp = -dpth*ones(1,NCells)-ones(Npings,1)*cellsize*double(0:1:NCells-1) - blockdis;
        elseif strcmp(direction,'up')
            dpth_temp = -dpth*ones(1,NCells)+ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis;
        end
        dpth_temp(find(m == 1)) = nan;
        dpth_temp_all = dpth_temp(:);
        
        % build time matrix
        time_temp = time*ones(1,NCells);
        time_temp(find(m==1)) = nan;
        time_temp_all = time_temp(:);
        
        % all velocity - range dependent/averaged
        velemcE_temp_all = velemc{1}(:);
        velemcN_temp_all = velemc{2}(:);
%         rangeaE_temp_all = rangea{1}(:);
%         rangeaN_temp_all = rangea{2}(:);
        
        %%
        z = (0:-boxs:-z_max)';
        
        temp_pro = nan(length(z),2);  % range dependdent velocity of each profile, first column: E/W; second column: N/S
%         temp_pro_ra = nan(length(z),2);  % range averaged velocity of each profile, first column: E/W; second column: N/S
        temp_time   = nan(length(z),1);  % box-averaged time
        for i = 1:length(z)
            index_box = [];
            index_box = find(z(i)-boxs/2<dpth_temp_all & dpth_temp_all<=z(i)+boxs/2);
            if length(index_box)>0
                temp_pro(i,1)  = nanmean(velemcE_temp_all(index_box));
                temp_pro(i,2)  = nanmean(velemcN_temp_all(index_box));
%                 temp_pro_ra(i,1)  = nanmean(rangeaE_temp_all(index_box));
%                 temp_pro_ra(i,2)  = nanmean(rangeaN_temp_all(index_box));
                temp_time(i,1)    = nanmean(time_temp_all(index_box));                       
            end
        end
        clear i
            
        vele_down = [vele_down temp_pro(:,1)];
        veln_down = [veln_down temp_pro(:,2)];
        time_down = [time_down temp_time];
        
        disp([filename,' total #= ',num2str(size(time_down,2))]);
    end
end

%% combine upcast and downcast
% time_up_av = nanmean(time_up,1);      % depth averaged time
% time_down_av = nanmean(time_down,1);  % depth averaged time
% 
% total_time_av = [time_up_av time_down_av];  % combine upcast with downcast
% total_time = [time_up time_down];
% total_vele = [vele_up vele_down];        
% total_veln = [veln_up veln_down];
% 
% [~, id] = sort(total_time_av);    % sort the time
% 
% out{1} = total_vele(:,id);
% out{2} = total_veln(:,id);
% out{3} = total_time(:,id);
% out{4} = -z'*ones(1,length(total_time_av));

end


