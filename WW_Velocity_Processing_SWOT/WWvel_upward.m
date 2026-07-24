function out = WWvel_upward(WWmeta,variables,upcast_bool,splitfiles,splitnum)
% out{1} = [];  % E-W velocity
% out{2} = [];  % N-s velocity
% out{3} = [];  % time
% out{4} = [];  % depth
% bofu zheng

%% set variables
if nargin<2
    num      = 1;   % 20 files in a group 
    blockdis = 0.1;  % blocking distance 0.1 m
    cellsize = 1; % cell size 1 m
    saprate  = 16;   % sampling rate 16 Hz
    boxs  = 0.5; % box size 0.5 m
    z_max    = 500;  % profile depth 500 m
    direction = 'up';
    disp('No varibles given, using default setup')
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
out{7} = [];  % Surface Vel
out{8} = [];  % Nav and temp
out{9} = [];  % Normalized amplitude
out{10} = [];  % Amplitude variance

dd0 = WWmeta.dd0;
if ((splitnum-1)*splitfiles+splitfiles)<ceil(length(WWmeta.dd0)/num)
    inds = (num*((splitnum-1)*splitfiles)+1):num:num*((splitnum-1)*splitfiles+splitfiles);  % starting index for each group
else
    inds = (num*((splitnum-1)*splitfiles)+1):num:num*(ceil(length(WWmeta.dd0)/num));  % starting index for each group
end

inde = inds(2:end)-1;

if (splitnum-1)*splitfiles+splitfiles<ceil(length(WWmeta.dd0)/num)
    inde = [inde num*((splitnum-1)*splitfiles+splitfiles)]; % ending index for each group
else
    inde = [inde length(WWmeta.dd0)]; % ending index for each group
end

%% analysis for upcast
Nav_Vars = {'Burst_WaterTemperature','Burst_Temperature','Burst_Heading','Burst_Pitch','Burst_Roll'}; % Specify nav variables to be bin averaged
for var = Nav_Vars;var=var{:};
    Nav.(var)=[];
end
vele_up = [];  % final product, E/W velocity
veln_up = [];  % final product, N/S velocity
velu_up = [];  % final product, Up/Down velocity
time_up = [];  % final product, time
shearu_up = []; % final product, E/W shear
shearv_up = []; % final product, N/S shear
surf_vel = []; % final product, Surface velocity
amp_up = []; % final product, normalized amplitude
amp_var_up = []; % final product, per bin amplitude std dev
num_per_bin = []; % num samples per bin
veln_up_var = [];  % N/S velocity variance
vele_up_var = [];  % E/W velocity variance
velu_up_var = [];  % Up/Down velocity variance
vele_corr_up = [];  % final product, E/W velocity, "Sail" corrected
veln_corr_up = [];  % final product, N/S velocity, "Sail" corrected

for q = 1:length(inds)
    filename = ['Profiles_upcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file
    load([WWmeta.propath_rearrange filename]);
    % find long enough profile
    index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_up))>0.2*z_max); % Change multiplier on z_max to change lengh of profile removed, multiplier represents fraction of max depth
     
    %% run through each upcast
    for k = 1:length(index)
        if isfield(AQDprofiles_up{index(k)},'Burst_Time')
            time=AQDprofiles_up{index(k)}.Burst_Time;
        elseif isfield(AQDprofiles_up{index(k)},'Burst_MatlabTimeStamp')
            time=AQDprofiles_up{index(k)}.Burst_MatlabTimeStamp;
        else
            disp("Error: Could not find a time variable")
        end
        
        dpth = AQDprofiles_up{index(k)}.Burst_Pressure;
        % Gather correlation data
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
        
        % Filter out low correlation samples
        m = corr1<50 | corr2<50 | corr3<50 | corr4<50;
        
        % Gather amplitude data
        amp=[];
        if isfield(AQDprofiles_up{index(k)},'Burst_AmpBeam1')
            amp(:,:,1) = double(AQDprofiles_up{index(k)}.Burst_AmpBeam1);
            amp(:,:,2) = double(AQDprofiles_up{index(k)}.Burst_AmpBeam2);
            amp(:,:,3) = double(AQDprofiles_up{index(k)}.Burst_AmpBeam3);
            amp(:,:,4) = double(AQDprofiles_up{index(k)}.Burst_AmpBeam4);
            amp = permute(amp,[1,3,2]);
        else
            amp = double(squeeze(AQDprofiles_up{index(k)}.Burst_Amplitude_Beam(:,:,:)));
        end
        
        % Gather velocity data
        beam_vel = [];
        if isfield(AQDprofiles_up{index(k)},'Burst_VelBeam1')
            beam_vel(:,:,1) = double(AQDprofiles_up{index(k)}.Burst_VelBeam1);
            beam_vel(:,:,2) = double(AQDprofiles_up{index(k)}.Burst_VelBeam2);
            beam_vel(:,:,3) = double(AQDprofiles_up{index(k)}.Burst_VelBeam3);
            beam_vel(:,:,4) = double(AQDprofiles_up{index(k)}.Burst_VelBeam4);
            beam_vel = permute(beam_vel,[1,3,2]);
        else
            beam_vel = AQDprofiles_up{index(k)}.Burst_Velocity_Beam;
        end
        
        % Sometimes there's a bug in aquisition, and amplitude/velocity
        % records are shorter than time records
        if size(amp,1)~=length(time) | size(beam_vel,1)~=length(time)
            continue
        end
        
        % Gather navigation data
        pitch = AQDprofiles_up{index(k)}.Burst_Pitch;
        roll = AQDprofiles_up{index(k)}.Burst_Roll;
        heading = AQDprofiles_up{index(k)}.Burst_Heading;
        if isfield(AQDprofiles_up{index(k)},'Burst_Accelerometer')
            acce = AQDprofiles_up{index(k)}.Burst_Accelerometer;
        elseif isfield(AQDprofiles_up{index(k)},'Burst_AccelerometerX') & ...
                isfield(AQDprofiles_up{index(k)},'Burst_AccelerometerY') & ...
                isfield(AQDprofiles_up{index(k)},'Burst_AccelerometerZ')
            acce = [AQDprofiles_up{index(k)}.Burst_AccelerometerX,...
                AQDprofiles_up{index(k)}.Burst_AccelerometerY,...
                AQDprofiles_up{index(k)}.Burst_AccelerometerZ];
        else
            disp('Error: could not find acceleration field')
        end
        
        if isfield(AQDprofiles_up{index(k)},'Burst_NCells')
            NCells = AQDprofiles_up{index(k)}.Burst_NCells(1); % number of cells
        elseif isfield(AQDprofiles_up{index(k)},'Burst_NumberofCells')
            NCells = AQDprofiles_up{index(k)}.Burst_NumberofCells(1);
        else
            disp('Error: could not find "Number of cells" field')
        end
        
        Npings = length(time);   % number of pings
        
        % Signuature 1000 beam geometry
        phi = [65,65,65,65]*pi/180;
        azi = [0,-90,180,90]*pi/180;
        
        % Get signature 1000 beam unit vectors, taking platform attitude
        % into account
        clear bX bY bZ
        for ibeam=1:4
          [bX(:,ibeam),bY(:,ibeam),bZ(:,ibeam)] = GetUnitVectors(phi(ibeam),azi(ibeam),roll*pi/180,pitch*pi/180);  
        end
        
        % caluate depth for each data point taking instrument orientaiton
        % and beam geometry into account
        ranges = (ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis + cellsize)/cosd(25);
        z_coords = zeros(size(ranges,1),size(ranges,2),4);
        for n = 1:size(ranges,1)
            if bZ(n,1)>0.1 & bZ(n,2)>0.1 & bZ(n,3)>0.1 & bZ(n,4)>0.1
                z_coords(n,:,:) = -dpth(n) + transpose(ranges(n,:)).*bZ(n,:);
            else
                z_coords(n,:,:) = nan(size(transpose(ranges(n,:)).*bZ(n,:)));
            end
        end
        
        % caluate the bin that best represents zero depth (surface echo)
        [min_dep,zero_bin] = min(abs(z_coords),[],2);
        zero_bin(min_dep>cellsize/2) = NaN;
        zero_bin = squeeze(zero_bin);
        
        dr = median(diff(ranges(1,:)));
        % Tranmission loss assuming circular spreading and 0.37dB/m
        % absorbtion
        TL = 10*log10((2*(ranges+dr/2)).^2)+2*0.37*(ranges+dr/2);
        % Normalize amplitude
        amp_norm = amp + permute(TL,[1,3,2]);
        
        amp_norm(:,:,1:2)=nan;
        amp_norm_1 = squeeze(amp_norm(:,1,:));
        amp_norm_2 = squeeze(amp_norm(:,2,:));
        amp_norm_3 = squeeze(amp_norm(:,3,:));
        amp_norm_4 = squeeze(amp_norm(:,4,:));
        
        % Use the beam geometry to calulate the expected sidelobe locations
        side_lobe_mask = zeros(size(beam_vel));
        for n = 1:size(beam_vel,1)
            [dist,side_lobe_range] = min(abs(dpth(n)-ranges(n,:)));
            if dist<1
                side_lobe_mask(n,:,side_lobe_range:end) = 1;
            end
        end
          
        % Isolate surface velocities for each beam (Velocity at zero depth)
        surface_beam_vel = nan(size(beam_vel,1),4);
        for b = 1:4
            for n = 1:size(beam_vel,1)
                if ~isnan(zero_bin(n,b))
                    surface_beam_vel(n,b) = beam_vel(n,b,zero_bin(n,b));
                end
            end
        end
        
        % remove data in the sidelobes
        beam_vel(smoothdata(side_lobe_mask,3,'movmean',3)>0)=NaN;
        
        % Preform rotation to put surface beam velocities in ENU
        % coordiantes
        ENU_surf = Beam2ENU(surface_beam_vel(:,1),surface_beam_vel(:,2),surface_beam_vel(:,3),surface_beam_vel(:,4), pitch, roll, heading);
        
        % build dpth matrix. This assumes a nominal instrument orientation.
        % We will sample beamwise velocities at these depths
        if strcmp(direction,'down')
            dpth_temp = -dpth*ones(1,NCells)-ones(Npings,1)*cellsize*double(0:1:NCells-1) - blockdis;
        elseif strcmp(direction,'up')
            dpth_temp = -dpth*ones(1,NCells)+ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis;
        end
        dpth_temp(find(m == 1)) = nan;
        % Reshape 2D depth arrays into 1D vectors
        dpth_temp_all = dpth_temp(:);
        dpth_temp_shear = dpth_temp(:,2:end-1);
        dpth_temp_shear_all = dpth_temp_shear(:);
        dpth_temp_all_1 = z_coords(:,:,1);dpth_temp_all_1=dpth_temp_all_1(:);
        dpth_temp_all_2 = z_coords(:,:,2);dpth_temp_all_2=dpth_temp_all_2(:);
        dpth_temp_all_3 = z_coords(:,:,3);dpth_temp_all_3=dpth_temp_all_3(:);
        dpth_temp_all_4 = z_coords(:,:,4);dpth_temp_all_4=dpth_temp_all_4(:);
        
        % Interpolate beam velocity data onto a even depth grid defined by the
        % nominal instrument orientation
        equal_depth_vel = nan(size(beam_vel));
        for b = 1:4
            for n = 1:size(beam_vel,1)
                if sum(~isnan(z_coords(n,:,b)))>1
                    equal_depth_vel(n,b,:) = interp1(z_coords(n,:,b),squeeze(beam_vel(n,b,:)),dpth_temp(n,:));
                end
            end
        end
        
        % Calulate velocity shear along each beam
        beamshear = permute((diff(permute(beam_vel(:,:,1:end-1),[1,3,2]),1,2)-fliplr(diff(fliplr(permute(beam_vel(:,:,2:end),[1,3,2])),1,2)))...
            ./(diff(z_coords(:,1:end-1,:),1,2)-fliplr(diff(fliplr(z_coords(:,2:end,:)),1,2))),[1,3,2]);
        
        % Interpolate beam shear data onto a even depth grid defined by the
        % nominal instrument orientation
        equal_depth_shear = nan(size(beam_vel,1),size(beam_vel,2),size(beam_vel,3)-2);
        for b = 1:4
            for n = 1:size(beam_vel,1)
                if sum(~isnan(z_coords(n,2:end-1,b)))>1
                    equal_depth_shear(n,b,:) = interp1(z_coords(n,2:end-1,b),squeeze(beamshear(n,b,:)),dpth_temp_shear(n,:));
                end
            end
        end
        
        MC_cell = WWcorr_beam(time,dpth,acce,pitch,roll,heading,equal_depth_vel,saprate)';
        equal_depth_vel_mc = MC_cell{1};
        
        ENU_surf(:,1) = ENU_surf(:,1) + MC_cell{3};
        ENU_surf(:,2) = ENU_surf(:,2) + MC_cell{4};
        ENU_surf(:,3) = ENU_surf(:,3) + MC_cell{5};
        
        % Get surface velocities
        surf = nanmean(ENU_surf);
        
        depth_steps = 1:max(ranges(:));
        % Attempt to predict wave orbital velocities from surface
        % velocities and remove them. I think it works alright, but can
        % cause bugs, so removed for now

        % wave_vel = zeros(size(ENU_surf,1),3,length(depth_steps));
        % for n = 1:3
        %     good_surf = find(~isnan(ENU_surf(:,n)));
        %     if isempty(good_surf)
        %         continue
        %     end
        % 
        %     if mod(good_surf(end)-good_surf(1)+1,2)==0
        %         good_surf(1) = good_surf(1)+1;
        %     end
        %     U_surf = fillmissing(ENU_surf(good_surf(1):good_surf(end),n),'linear');
        %     f_surf = saprate*(0:(length(U_surf)/2))/length(U_surf);
        %     f_surf = [f_surf,-fliplr(f_surf(2:end))];
        %     k_surf = (2*pi*f_surf).^2/9.8;
        % 
        %     dep_atten = exp(-transpose(k_surf).*depth_steps);
        % 
        %     U_surf = detrend(U_surf,'omitnan');
        %     U_surf_f = [flip(U_surf);U_surf;flip(U_surf)];
        %     [b,a] = butter(7,[0.07 1.2]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
        %     U_surf_filt = filtfilt(b,a,double(U_surf_f));
        %     U_surf_filt = U_surf_filt(length(U_surf_filt)*1/3+1:length(U_surf_filt)*2/3);  % chunk back to original length
        % 
        %     ft_surf = fft(U_surf_filt);
        %     surf_atten = real(ifft(ft_surf.*dep_atten));
        % 
        %     wave_vel(good_surf(1):good_surf(end),n,:) = permute(surf_atten,[1,3,2]);
        % end
        % 
        % wave_vel_interp = zeros(size(dpth_temp,1),3,size(dpth_temp,2));
        % [xs,zs] = meshgrid(1:length(wave_vel),-depth_steps);
        % for n = 1:3
        %     wave_vel_interp(:,n,:) = permute(interp2(xs,zs,squeeze(wave_vel(:,n,:))',repmat((1:size(dpth_temp,1))',1,size(dpth_temp,2)),double(dpth_temp)),[1,3,2]);
        % end
            
        % Rotate shear and velocity from beam coordianates to ENU
        % coordinates
        equal_depth_shear_ENU = Beam2ENU(squeeze(equal_depth_shear(:,1,:)),squeeze(equal_depth_shear(:,2,:)),...
            squeeze(equal_depth_shear(:,3,:)),squeeze(equal_depth_shear(:,4,:)), pitch, roll, heading);
        
        equal_depth_ENU = Beam2ENU(squeeze(equal_depth_vel_mc(:,1,:)),squeeze(equal_depth_vel_mc(:,2,:)),...
            squeeze(equal_depth_vel_mc(:,3,:)),squeeze(equal_depth_vel_mc(:,4,:)), pitch, roll, heading);
        
        % To do wave orbital removal this must be un-commented
        % wave_vel_interp(isnan(wave_vel_interp))=0;
%         vele{1} = squeeze(equal_depth_ENU(:,1,:))-squeeze(wave_vel_interp(:,1,:));
%         vele{2} = squeeze(equal_depth_ENU(:,2,:))-squeeze(wave_vel_interp(:,2,:));
%         vele{3} = squeeze(equal_depth_ENU(:,3,:))-squeeze(wave_vel_interp(:,3,:));
        
        % Correct for horizontal motion of the wirewalker on the wire using
        % the veritcal speeed and orientation of the wirewalker
        if variables.sail_corr==1
            dp_dt = -diff(dpth)./(diff(time)*86400);
            dp_dt_sm = smoothdata(dp_dt,'gaussian',24);
            dp_dt_sm = [0;dp_dt_sm];

            % Get unit-vector for the +z axis on the wirewalker in ENU coodinates
            ENU = XYZ2ENU(ones(size(pitch))*variables.z_unit(1),ones(size(pitch))*variables.z_unit(2),ones(size(pitch))*variables.z_unit(3),pitch,roll,heading);
        
            % Horizontal unit vector
            b_h = sqrt(ENU(:,1).^2+ENU(:,2).^2);
            
            b_x = ENU(:,1);
            b_y = ENU(:,2);
            b_z = ENU(:,3);
            
            % Horizontal velocity compnent from wirewalker rise rate
            v_h = dp_dt_sm.*b_h./b_z;           
            % N-S velocity component
            v_x = v_h.*b_x./b_h;         
            % E-W velocity component
            v_y = v_h.*b_y./b_h;
        else
            v_x = zeros(size(time));
            v_y = zeros(size(time));
        end

        vele{1} = squeeze(equal_depth_ENU(:,1,:));
        vele{2} = squeeze(equal_depth_ENU(:,2,:));
        vele{3} = squeeze(equal_depth_ENU(:,3,:));
        
        vele_corr{1} = squeeze(equal_depth_ENU(:,1,:))+v_x;
        vele_corr{2} = squeeze(equal_depth_ENU(:,2,:))+v_y;
        
        % Smooth shears in time covering ~1-2 secs of data (0.5-1.5m equivelent vertical smoothing)
        shearu1 = smoothdata(smoothdata(squeeze(equal_depth_shear_ENU(:,1,:)),2,'movmean',1),'gaussian',8);
        shearv1 = smoothdata(smoothdata(squeeze(equal_depth_shear_ENU(:,2,:)),2,'movmean',1),'gaussian',8);
        % Shears near sonar head are bad (stagnation point?)
        shearu1(:,1:2)= NaN;
        shearv1(:,1:2)= NaN;

        % remove velocity that has low correlation
        vele{1}(find(m == 1)) = nan;
        vele{2}(find(m == 1)) = nan;
        vele{3}(find(m == 1)) = nan;
        
        vele{3}(repmat(smoothdata(MC_cell{5},'gaussian',30)<0.8*median(MC_cell{5}),1,size(vele{3},2))) = NaN;
        
        % build time matrix
        time_temp = time*ones(1,NCells);
        time_temp(find(m==1)) = nan;
        time_temp_all = time_temp(:);

        % reshape velocity and shear from 2D arrays into vectors
        velemcE_temp_all = vele{1}(:);  % E/W
        velemcN_temp_all = vele{2}(:);  % N/S
        velemcE_corr_temp_all = vele_corr{1}(:);  % E/W
        velemcN_corr_temp_all = vele_corr{2}(:);  % N/S
        velemcU_temp_all = vele{3}(:);  % N/S
        shearu1_temp = shearu1(:);
        shearv1_temp = shearv1(:);
        
        %% box averaging
        z = (0:-boxs:-z_max)';
        % Temporary arrays to store binned data
        temp_pro = nan(length(z),3);  % each profile, first column: E/W; second column: N/S
        temp_pro_corr = nan(length(z),2);  % each profile, first column: E/W; second column: N/S
        temp_pro_amp = nan(length(z),4); % each profile, 1 column per beam
        temp_pro_amp_var = nan(length(z),4); % each profile, 1 column per beam
        temp_pro_shear = nan(length(z),2); % each profile, first column: E/W; second column: N/S
        temp_time = nan(length(z),1);  % box-averaged time        
        temp_bin_num = nan(length(z),1); %number of samples per bin
        temp_pro_vel_var = nan(length(z),3);  % each profile variance per bin, first column: E/W; second column: N/S

        for var = Nav_Vars;var=var{:};
            temp_pro_Nav.(var) = nan(length(z),1); % each profile, 1 field per nav variable
        end
        
        for i = 1:length(z)
            % get samples within a given depth bin for nominal instrument orientation 
            index_box = find(z(i)-boxs/2<dpth_temp_all & dpth_temp_all<=z(i)+boxs/2);
            % same, but with beamwise depths
            index_box_1 = find(z(i)-boxs/2<dpth_temp_all_1 & dpth_temp_all_1<=z(i)+boxs/2);
            index_box_2 = find(z(i)-boxs/2<dpth_temp_all_2 & dpth_temp_all_2<=z(i)+boxs/2);
            index_box_3 = find(z(i)-boxs/2<dpth_temp_all_3 & dpth_temp_all_3<=z(i)+boxs/2);
            index_box_4 = find(z(i)-boxs/2<dpth_temp_all_4 & dpth_temp_all_4<=z(i)+boxs/2);
            % and same, but for 1D variables (Nav)
            index_box_1D = find(z(i)-boxs/2<-dpth & -dpth<=z(i)+boxs/2);
            
            % Bin amplitudes, and amplitude variance
            if ~isempty(index_box_1) & ~isempty(index_box_2) && ~isempty(index_box_3) & ~isempty(index_box_4)
                temp_pro_amp(i,1)  = nanmean(amp_norm_1(index_box_1));
                temp_pro_amp(i,2)  = nanmean(amp_norm_2(index_box_2));
                temp_pro_amp(i,3)  = nanmean(amp_norm_3(index_box_3));
                temp_pro_amp(i,4)  = nanmean(amp_norm_4(index_box_4)); 
                
                temp_pro_amp_var(i,1)  = nanstd(amp_norm_1(index_box_1));
                temp_pro_amp_var(i,2)  = nanstd(amp_norm_2(index_box_2));
                temp_pro_amp_var(i,3)  = nanstd(amp_norm_3(index_box_3));
                temp_pro_amp_var(i,4)  = nanstd(amp_norm_4(index_box_4)); 
            end
            
            % Bin time and velocity
            if length(index_box)>0 & sum(~isnan(velemcE_temp_all(index_box)))>10 &...
                    sum(~isnan(velemcN_temp_all(index_box)))>10
                temp_pro(i,1)  = nanmean(velemcE_temp_all(index_box));
                temp_pro(i,2)  = nanmean(velemcN_temp_all(index_box));
                temp_pro(i,3)  = nanmean(velemcU_temp_all(index_box));  
                temp_pro_corr(i,1)  = nanmean(velemcE_corr_temp_all(index_box));
                temp_pro_corr(i,2)  = nanmean(velemcN_corr_temp_all(index_box));
                temp_time(i,1) = nanmean(time_temp_all(index_box));
                temp_bin_num(i,1)  = sum(~isnan(velemcE_temp_all(index_box)) & ~isnan(velemcN_temp_all(index_box)));
                temp_pro_vel_var(i,1)  = nanstd(velemcE_temp_all(index_box));
                temp_pro_vel_var(i,2)  = nanstd(velemcN_temp_all(index_box)); 
                temp_pro_vel_var(i,3)  = nanstd(velemcU_temp_all(index_box));      
            end
            % Bin nav variables
            if length(index_box_1D)>0
                for var = Nav_Vars;var=var{:};
                    try
                        temp_pro_Nav.(var)(i) = nanmean(AQDprofiles_up{index(k)}.(var)(index_box_1D));
                    catch
                        continue
                    end
                end
            end
        end
        clear i
        
        % Bin shear separately (Different z coordinates) 
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
        
        % Add profile to profiles already calulated.
        vele_up = [vele_up temp_pro(:,1)];
        veln_up = [veln_up temp_pro(:,2)];
        velu_up = [velu_up temp_pro(:,3)];
        vele_corr_up = [vele_corr_up temp_pro_corr(:,1)];
        veln_corr_up = [veln_corr_up temp_pro_corr(:,2)];
        amp_up = cat(3,amp_up,temp_pro_amp);
        amp_var_up = cat(3,amp_var_up,temp_pro_amp_var);
        shearu_up = [shearu_up temp_pro_shear(:,1)];
        shearv_up = [shearv_up temp_pro_shear(:,2)];
        time_up = [time_up temp_time];
        surf_vel = [surf_vel surf(1:2)'];
        num_per_bin = [num_per_bin temp_bin_num];
        vele_up_var = [vele_up_var temp_pro_vel_var(:,1)];
        veln_up_var = [veln_up_var temp_pro_vel_var(:,2)];
        velu_up_var = [velu_up_var temp_pro_vel_var(:,3)];

        for var = Nav_Vars;var=var{:};
            Nav.(var) = [Nav.(var),temp_pro_Nav.(var)];
        end
        
        disp(['Processing Mean Velocity ' filename,' total #= ',num2str(size(time_up,2))]);
    end
end

%% Add to outputs
out{1} = vele_up;
out{2} = veln_up;
out{3} = time_up;
out{4} = -z*ones(1,size(time_up,2));
out{5} = shearu_up;
out{6} = shearv_up;
out{7} = surf_vel;
out{8} = Nav;
out{9} = permute(amp_up,[1,3,2]);
out{10} = permute(amp_var_up,[1,3,2]);
out{11} = num_per_bin;
out{12} = vele_up_var;
out{13} = veln_up_var;
out{14} = velu_up;
out{15} = velu_up_var;
out{16} = vele_corr_up;
out{17} = veln_corr_up;

%% analysis for downcast - still needs work to be done
if upcast_bool == 1

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
time_up_av = nanmean(time_up,1);      % depth averaged time
time_down_av = nanmean(time_down,1);  % depth averaged time

total_time_av = [time_up_av time_down_av];  % combine upcast with downcast
total_time = [time_up time_down];
total_vele = [vele_up vele_down];        
total_veln = [veln_up veln_down];

[~, id] = sort(total_time_av);    % sort the time

out{1} = total_vele(:,id);
out{2} = total_veln(:,id);
out{3} = total_time(:,id);
out{4} = -z'*ones(1,length(total_time_av));

end


