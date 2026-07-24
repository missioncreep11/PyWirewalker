max_depth = 110;      % Depth of deepest ADCP bin for whole deployment
min_depth = 0;        % Depth of shallowest ADCP bin for whole deployment
time_bin = 60;        % Size of final time bins in seconds
box_size = 1;         % Size of final depth bins in meters
cellsize = 1;         % Size of cell in ADCP setup
blockdis = 0.5;       % blanking distance from ADCP setup
saprate = 4;          % Sample rate (Hz) for ADCP setup
direction = 'up';     % Look direction. Only "up" has been tested
dep_start = datenum([2025,3,31,21,38,0]);   % Deployment period start
dep_end = datenum([2025,4,16,0,0,0]);       % Deployment period end
raw_data_files = dir('/Users/Devon/Documents/GradSchool/DTS/TLC25_Sig500/Raw/*.mat'); % Path to raw files
savefile = '/Users/Devon/Documents/GradSchool/DTS/TLC25_Sig500/Binned_ADCP.mat';      % path to save file

% initalize binned outputs
velE_binned = zeros(length(min_depth:box_size:max_depth),length(dep_start:time_bin/86400:dep_end),2);
velN_binned = zeros(length(min_depth:box_size:max_depth),length(dep_start:time_bin/86400:dep_end),2);
velU_binned = zeros(length(min_depth:box_size:max_depth),length(dep_start:time_bin/86400:dep_end),2);
amp_binned = zeros(length(min_depth:box_size:max_depth),length(dep_start:time_bin/86400:dep_end),4,2);

for i = 28:length(raw_data_files)
        % Load data
        load([raw_data_files(i).folder '/' raw_data_files(i).name]);        
        disp(['Processing file ' num2str(i)])
        
        % Gather time data
        if isfield(Data,'Burst_Time')
            time=Data.Burst_Time;
        elseif isfield(Data,'Burst_MatlabTimeStamp')
            time=Data.Burst_MatlabTimeStamp;
        else
            disp("Error: Could not find a time variable")
        end
        
        % Gather pressure data
        dpth = Data.Burst_Pressure;

        % Gather correlation data
        if isfield(Data,'Burst_CorBeam1')
            corr1 = double(Data.Burst_CorBeam1);
            corr2 = double(Data.Burst_CorBeam2);
            corr3 = double(Data.Burst_CorBeam3);
            corr4 = double(Data.Burst_CorBeam4);
        else
            corr1 = double(squeeze(Data.Burst_Correlation_Beam(:,1,:)));
            corr2 = double(squeeze(Data.Burst_Correlation_Beam(:,2,:)));
            corr3 = double(squeeze(Data.Burst_Correlation_Beam(:,3,:)));
            corr4 = double(squeeze(Data.Burst_Correlation_Beam(:,4,:)));
        end
        
        % Filter out low correlation samples
        m = corr1<50 | corr2<50 | corr3<50 | corr4<50;
        
        % Gather amplitude data
        amp=[];
        if isfield(Data,'Burst_AmpBeam1')
            amp(:,:,1) = double(Data.Burst_AmpBeam1);
            amp(:,:,2) = double(Data.Burst_AmpBeam2);
            amp(:,:,3) = double(Data.Burst_AmpBeam3);
            amp(:,:,4) = double(Data.Burst_AmpBeam4);
            amp = permute(amp,[1,3,2]);
        else
            amp = double(squeeze(Data.Burst_Amplitude_Beam(:,:,:)));
        end
        
        % Gather velocity data
        beam_vel = [];
        if isfield(Data,'Burst_VelBeam1')
            beam_vel(:,:,1) = double(Data.Burst_VelBeam1);
            beam_vel(:,:,2) = double(Data.Burst_VelBeam2);
            beam_vel(:,:,3) = double(Data.Burst_VelBeam3);
            beam_vel(:,:,4) = double(Data.Burst_VelBeam4);
            beam_vel = permute(beam_vel,[1,3,2]);
        else
            beam_vel = Data.Burst_Velocity_Beam;
        end
        
        % Gather navigation data
        pitch = Data.Burst_Pitch;
        roll = Data.Burst_Roll;
        heading = Data.Burst_Heading;
        if isfield(Data,'Burst_Accelerometer')
            acce = Data.Burst_Accelerometer;
        elseif isfield(Data,'Burst_AccelerometerX') & ...
                isfield(Data,'Burst_AccelerometerY') & ...
                isfield(Data,'Burst_AccelerometerZ')
            acce = [Data.Burst_AccelerometerX,...
                Data.Burst_AccelerometerY,...
                Data.Burst_AccelerometerZ];
        else
            disp('Error: could not find acceleration field')
        end
        
        if isfield(Data,'Burst_NCells')
            NCells = Data.Burst_NCells(1); % number of cells
        elseif isfield(Data,'Burst_NumberofCells')
            NCells = Data.Burst_NumberofCells(1);
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
        
        dr = median(diff(ranges(1,:)));
        % Tranmission loss assuming circular spreading and 0.13dB/m
        % absorbtion
        TL = 10*log10((2*(ranges+dr/2)).^2)+2*0.13*(ranges+dr/2);
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
          
        % remove data in the sidelobes
        beam_vel(smoothdata(side_lobe_mask,3,'movmean',3)>0)=NaN;
        
        % build dpth matrix. This assumes a nominal instrument orientation.
        % We will sample beamwise velocities at these depths       
        dpth_temp = -dpth*ones(1,NCells)+ones(Npings,1)*cellsize*double(0:1:NCells-1) + blockdis;
        dpth_temp(find(m == 1)) = nan;
  
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
        
        % Motion correct beam data
        MC_cell = WWcorr_beam(time,dpth,acce,pitch,roll,heading,equal_depth_vel,saprate)';
        equal_depth_vel_mc = MC_cell{1};     
        
        % rotate into beam coordinates
        equal_depth_ENU = Beam2ENU(squeeze(equal_depth_vel_mc(:,1,:)),squeeze(equal_depth_vel_mc(:,2,:)),...
            squeeze(equal_depth_vel_mc(:,3,:)),squeeze(equal_depth_vel_mc(:,4,:)), pitch, roll, heading);
        
        vele{1} = squeeze(equal_depth_ENU(:,1,:));
        vele{2} = squeeze(equal_depth_ENU(:,2,:));
        vele{3} = squeeze(equal_depth_ENU(:,3,:));

        % remove velocity that has low correlation
        vele{1}(find(m == 1)) = nan;
        vele{2}(find(m == 1)) = nan;
        vele{3}(find(m == 1)) = nan;
        
        %%
        
        % smooth over 20 secs to remove wave signal
        vele{1} = smoothdata(vele{1},'gaussian',saprate*20);
        vele{2} = smoothdata(vele{2},'gaussian',saprate*20);
        vele{3} = smoothdata(vele{3},'gaussian',saprate*20);
        
        % Bin velocity and amplitude data
        depth_ind = ceil((dpth_temp+max_depth)/box_size);
        time_ind = ceil((time - dep_start)/(time_bin/86400));
        time_ind = repmat(time_ind,1,size(depth_ind,2));

        for n = 1:numel(depth_ind)
            if ~isnan(depth_ind(n)) & ~isnan(time_ind(n)) & ...
                depth_ind(n)>0 & time_ind(n)>0 & ...            
                depth_ind(n)<size(velE_binned,1) & time_ind(n)<size(velE_binned,2) & ...
                ~isnan(vele{1}(n)) & ~isnan(vele{2}(n)) & ~isnan(vele{3}(n))

                velN_binned(depth_ind(n),time_ind(n),1) = velN_binned(depth_ind(n),time_ind(n),1) + vele{2}(n);
                velN_binned(depth_ind(n),time_ind(n),2) = velN_binned(depth_ind(n),time_ind(n),2) + 1;
    
                velE_binned(depth_ind(n),time_ind(n),1) = velE_binned(depth_ind(n),time_ind(n),1) + vele{1}(n);
                velE_binned(depth_ind(n),time_ind(n),2) = velE_binned(depth_ind(n),time_ind(n),2) + 1;
    
                velU_binned(depth_ind(n),time_ind(n),1) = velU_binned(depth_ind(n),time_ind(n),1) + vele{3}(n);
                velU_binned(depth_ind(n),time_ind(n),2) = velU_binned(depth_ind(n),time_ind(n),2) + 1;
            end

            if ~isnan(depth_ind(n)) & ~isnan(time_ind(n)) & ...
                depth_ind(n)>0 & time_ind(n)>0 & ...            
                depth_ind(n)<size(velE_binned,1) & time_ind(n)<size(velE_binned,2) & ...
                ~isnan(amp_norm_1(n)) & ~isnan(amp_norm_2(n)) & ~isnan(amp_norm_3(n)) & ~isnan(amp_norm_4(n))

                amp_binned(depth_ind(n),time_ind(n),1,1) = amp_binned(depth_ind(n),time_ind(n),1,1) + amp_norm_1(n);
                amp_binned(depth_ind(n),time_ind(n),1,2) = amp_binned(depth_ind(n),time_ind(n),1,2) + 1;
    
                amp_binned(depth_ind(n),time_ind(n),2,1) = amp_binned(depth_ind(n),time_ind(n),2,1) + amp_norm_2(n);
                amp_binned(depth_ind(n),time_ind(n),2,2) = amp_binned(depth_ind(n),time_ind(n),2,2) + 1;
    
                amp_binned(depth_ind(n),time_ind(n),3,1) = amp_binned(depth_ind(n),time_ind(n),3,1) + amp_norm_3(n);
                amp_binned(depth_ind(n),time_ind(n),3,2) = amp_binned(depth_ind(n),time_ind(n),3,2) + 1;

                amp_binned(depth_ind(n),time_ind(n),4,1) = amp_binned(depth_ind(n),time_ind(n),4,1) + amp_norm_4(n);
                amp_binned(depth_ind(n),time_ind(n),4,2) = amp_binned(depth_ind(n),time_ind(n),4,2) + 1;
            end
        end
end

% average velocity and amplitude data
velN_binned = velN_binned(:,:,1)./velN_binned(:,:,2);
velE_binned = velE_binned(:,:,1)./velE_binned(:,:,2);
velU_binned = velU_binned(:,:,1)./velU_binned(:,:,2);

amp_binned = amp_binned(:,:,:,1)./amp_binned(:,:,:,2);

% build and save output structure
ADCP.depth = transpose((max_depth:-box_size:min_depth) - box_size/2);
ADCP.time = (dep_start:time_bin/86400:dep_end) + time_bin/86400/2;
ADCP.velE = velE_binned;
ADCP.velN = velN_binned;
ADCP.velU = velU_binned;
ADCP.amp = amp_binned;

save(savefile,'ADCP')
