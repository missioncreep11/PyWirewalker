function turb = WWturb_upward(WWmeta,variables,splitfiles,splitnum)

%% set variables
num      = variables.NUM_combining_files;
blockdis = variables.HRblockdis;
cellsize = variables.HRcellsize;
boxs  = variables.HRboxsize;
z_max    = variables.z_max;
HRbeams = variables.HRbeams;
w = warning ('off','all');

%%

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
turb.ep = [];turb.N = [];turb.SNR = [];ind=1;turb.A = [];
turb.ep_struct = [];turb.N_struct = [];turb.time = [];turb.A_struct = [];
turb.spec = [];turb.struct_fun = [];turb.r = [];turb.k = [];
turb.corr = [];turb.spec_num = [];turb.slope=[];turb.N_slope=[];
for q = 1:length(inds)
    filename = ['Profiles_upcast_',WWmeta.name_aqd,'_',num2str(inds(q)),'_',num2str(inde(q)),'.mat']; % Name of the processed file
    load([WWmeta.propath_rearrange filename]);
    % find long enough profile
    index=find(abs(cellfun(@(x) x.Burst_Pressure(end)-x.Burst_Pressure(1),AQDprofiles_up))>0.2*z_max); % Change multiplier on z_max to change lengh of profile removed, multiplier represents fraction of max depth
     
    %% run through each upcast
    for k = 1:length(index)
        profile = AQDprofiles_up{index(k)};
        
        if ~isfield(profile,'IBurstHR_Pressure')
            continue
        end
        
        % Get depth and vertical motion of the platform
        depth = profile.IBurstHR_Pressure;
        if isfield(profile,'IBurstHR_MatlabTimeStamp')
            time = profile.IBurstHR_MatlabTimeStamp;
        elseif isfield(profile,'IBurstHR_Time')
            time = profile.IBurstHR_Time;
        end
        if isfield(profile,'IBurstHR_Ambiguity')
            v_a = profile.IBurstHR_Ambiguity;
        elseif isfield(profile,'IBurstHR_AmbiguityVel')
            v_a = profile.IBurstHR_AmbiguityVel;
        end
        if isfield(profile,'IBurstHR_NumberofCells')
            N_Cells = profile.IBurstHR_NumberofCells;
        elseif isfield(profile,'IBurstHR_NCells')
            N_Cells = profile.IBurstHR_NCells;
        end
        
        roll = profile.IBurstHR_Roll;
        pitch = profile.IBurstHR_Pitch;
                
        % disp(['Processifng file ' num2str(q) ' profile ' num2str(index(k)) ' ' datestr(mean(time,'omitnan))])
        disp(['Processing HR turbulence ' filename,' total #= ',num2str(ind)]);

        % Signuature 1000 beam geometry
        phi = [65,65,65,65,90]*pi/180;
        azi = [0,-90,180,90,0]*pi/180;
        
        struct_res = boxs; % true resolution, tubulence will be sampled twice per resolution interval
        Afunc = @(rin) [rin./rin,rin.^(2/3)]; % Struture function formulat to fit
        vel_inds = round(struct_res/cellsize/2); % length of structure functions = resolution/2
            
        z_grid = 0:cellsize:z_max;
        bin_centers = vel_inds*cellsize/2:vel_inds*cellsize:z_max;
        
        for iBeam = 1:length(HRbeams)
            % velocity and correltion
            vel = eval(['profile.IBurstHR_VelBeam' num2str(HRbeams(iBeam))]);
            corr = eval(['profile.IBurstHR_CorBeam' num2str(HRbeams(iBeam))]);

            if size(vel,1)~=length(depth)
                disp('Mismatch between velocity and pressure data')
                continue
            end
    
            % Remove low correlation velocity values
            vel_masked = vel;
            vel_masked(corr<50)=NaN;
    
            % Cast all velocities as angles, remove the mean. This elimiates most
            % large spikes, which come from incorrect alias corections
            vel_angle = (vel_masked)/(median(v_a,'omitnan'));
            vel_angle = vel_angle*2*pi;
            angle_zeroed_cmplx = exp(1i*vel_angle);
            angle_zeroed_cmplx = angle_zeroed_cmplx.*exp(-1i*angle(mean(angle_zeroed_cmplx(:,15:end-1),2,'omitnan')));
            
            % Calulate the vertical position of each bin
            [bX,bY,bZ] = GetUnitVectors(phi(HRbeams(iBeam)),azi(HRbeams(iBeam)),roll*pi/180,pitch*pi/180);
    
            z = blockdis:cellsize:cellsize*double(median(N_Cells))+blockdis-cellsize;
            range = blockdis:cellsize:cellsize*double(median(N_Cells))+blockdis-cellsize;
            z_coord = depth-transpose(bZ).*repmat(z,length(depth),1);
    
            % demean velocity. First of sevral demeaning steps
            vel_zeroed = angle(angle_zeroed_cmplx)/(2*pi)*median(v_a,'omitnan');
    
            % find point where the stagnation point "ends". This is defined as the
            % point where the median velocity over a profile gets within 4mm/s of 
            % the median value along that profile. Any velocity before this point
            % will be removed. Any further stagnation point effects will be removed
            % with a linear detrend.
            avg_deep_prof = median(vel_zeroed(depth>20,:),'omitnan');
            cutoff_vel = median(avg_deep_prof,'omitnan')+0.004;
            cutoff_ind = find(diff(avg_deep_prof>cutoff_vel)==-1);
            if isempty(cutoff_ind)
                continue
            end
            cutoff_ind = cutoff_ind(end);
    
            % Remove velocity "spikes", single points with large dv/dz before and
            % after
            vel_zeroed_filled = fillmissing(vel_zeroed,'constant',0);
            down_spike = (int16(diff(vel_zeroed_filled(:,1:end-1),1,2)<-0.005) + int16(diff(vel_zeroed_filled(:,2:end),1,2)>0.005))==2;
            up_spike = (int16(diff(vel_zeroed_filled(:,1:end-1),1,2)>0.005) + int16(diff(vel_zeroed_filled(:,2:end),1,2)<-0.005))==2;
            down_spike = [zeros(size(down_spike,1),1),down_spike];
            up_spike = [zeros(size(up_spike,1),1),up_spike];
            vel_zeroed(down_spike==1)=NaN;
            vel_zeroed(up_spike==1)=NaN;
    
            % de-mean velocities
            vel_zeroed_use = transpose(detrend(transpose(vel_zeroed(:,cutoff_ind:end-2)),'omitnan'));
            z_coord_use = z_coord(:,cutoff_ind:end-2);
            corr_use =  corr(:,cutoff_ind:end-2);
    
            % Structure Function Method
            range_sub = cellsize:cellsize:cellsize*vel_inds; % range axis for structure function fits
            D_save = nan(size(z_grid,2),vel_inds);
            for g = 1:size(z_grid,2)
                % find closest velocity estimate to grid level
                [min_val,i] = min(abs(z_coord_use-z_grid(g)),[],2);
                % find all velocity profiles containing a grid point
                dep_mask = min_val<=cellsize/2;
                i = i(dep_mask);
                % subset velocity
                vel_sub = vel_zeroed_use(dep_mask,:);
    
                D_tot = zeros(2,vel_inds);
                for j = 1:sum(dep_mask)
                    if i(j)>1 & i(j)<size(vel_sub,2)
                        % First side of structure function
                        D1 = (vel_sub(j,i(j)) - vel_sub(j,i(j)+1:min([i(j)+vel_inds,size(vel_sub,2)]))).^2;
                        % second side of structure function
                        D2 = fliplr((vel_sub(j,i(j)) - vel_sub(j,max([i(j)-vel_inds,1]):i(j)-1)).^2);
                        nanmask_D1 = isnan(D1);
                        nanmask_D2 = isnan(D2);
                        D1(nanmask_D1) = 0;
                        D2(nanmask_D2) = 0;
                        % Average structure functions
                        D_tot(1,1:length(D1)) = D_tot(1,1:length(D1)) + D1;
                        D_tot(1,1:length(D2)) = D_tot(1,1:length(D2)) + D2;
                        D_tot(2,1:length(D1)) = D_tot(2,1:length(D1)) + ~nanmask_D1;
                        D_tot(2,1:length(D2)) = D_tot(2,1:length(D2)) + ~nanmask_D2;
                    end
                end
    
                D_save(g,:) = D_tot(1,:)./D_tot(2,:);
            end
    
            % Smooth structure function over one resolution bin
            D_save_sm = smoothdata(D_save,'gaussian',vel_inds*2);
            % Decimate structure function every hald resolution bin
            D_save_dec = D_save_sm(floor(vel_inds/2+1):vel_inds:end,:);
            z_grid_dec = z_grid(floor(vel_inds/2+1):vel_inds:end);
    
            if size(D_save_dec,2)>10
                fit_end = 15;
            else
                fit_end = size(D_save_dec,2);
            end
    
            % fit structure functions
            A_turb = nan(size(D_save_dec,1),1);
            N = nan(size(D_save_dec,1),1);
            A = Afunc(transpose(range_sub(1:fit_end)));
            for i=1:size(D_save_dec,1)
                Coeffs = (A'*A)^-1*A'*transpose(D_save_dec(i,1:fit_end));
                N(i) = Coeffs(1);A_turb(i) = Coeffs(2);
            end
    
            % calculate epsilon
            ep = real((A_turb/2).^(3/2));
            
            turb.ep_struct(:,ind,iBeam) = transpose(ep);
            turb.A_struct(:,ind,iBeam) = transpose(A_turb);
            turb.N_struct(:,ind,iBeam) = transpose(N);
            turb.struct_fun(:,ind,:,iBeam) = permute(D_save_dec,[1,3,2]);
    
            % Wavenumber spectrum method
            ks_len = 100;
            ks = 1/cellsize*(0:ks_len/2-1)/ks_len;
            dk = 1/cellsize/ks_len;
            k_bin_edges = [(ks(2:end)+ks(1:end-1))/2 - dk;(ks(2:end)+ks(1:end-1))/2];
            k_bin_edges(1,end+1) = k_bin_edges(2,end);
            k_bin_edges(2,end) = k_bin_edges(1,end)+dk;
    
            for i = 1:length(bin_centers)
                dep_mask = z_coord_use>=bin_centers(i)-boxs/2 & z_coord_use<bin_centers(i)+boxs/2;
                corr_avg = mean(corr_use(dep_mask),'omitnan');
                
                prof_num = find(sum(dep_mask,2)>4);
                spec = nan(length(ks),length(prof_num));
                ind_spec = 1;
                for j = prof_num'
                    vel_subset = vel_zeroed_use(j,dep_mask(j,:));
                    vel_subset = fillmissing(vel_subset-mean(vel_subset,'omitnan'),'constant',0);
    
                    window = hamming(length(vel_subset));
                    ft = abs(fft(vel_subset'.*window)).^2*(2/sum(window.^2))/length(vel_subset);
                    f = 1/cellsize*(0:(length(vel_subset)/2))/length(vel_subset);
                    df = median(2*pi*diff(f),'omitnan');
                    ft = ft/df;
    
                    spec(:,ind_spec) = interp1(f(2:end),transpose(ft(2:length(f))),ks);
                    ind_spec = ind_spec+1;
                end
    
                num_specs = max(sum(~isnan(spec),2));
                bad_specs = sum(isnan(spec))==size(spec,1) | sum(spec,'omitnan')>mean(sum(spec,'omitnan'),'omitnan')+2*std(sum(spec,'omitnan'),'omitnan');
                spec(:,bad_specs)=[];
    
                spec_fit = mean(spec,2,'omitnan');
                avg_num = sum(~isnan(spec),2,'omitnan');
                spec_fit(spec_fit<10^-8) = NaN;
                spec_fit(avg_num<max(avg_num)/2.5)=NaN;
        
                fit_cut = find(ks<0.5,1,'last')+1;
                A = FitKolmogorov(2*pi*ks(fit_cut:end),spec_fit(fit_cut:end));

                 spec_dN = real(log10(spec_fit-A(1)));
                spec_dN(abs(imag(log10(spec_fit-A(1))))>0) = NaN;
                numel_spec = length(spec_fit)-fit_cut+1;
                cutoff = find(spec_fit(2:end)<1.1*A(1),1);
                if isempty(cutoff)
                    cutoff = length(spec_fit);
                end

                SLfunc = @(kin) [kin./kin,kin];
                k_use = log10(2*pi*ks(fit_cut:cutoff));spec_use = spec_dN(fit_cut:cutoff);
                nanmask = isnan(spec_use);
                k_use(nanmask)=[];spec_use(nanmask)=[];
                if isempty(k_use)
                    Coeffs_SL = [NaN,NaN];
                else
                    A_SL = SLfunc(k_use');
                    Coeffs_SL = (A_SL'*A_SL)^-1*A_SL'*spec_use;
                end
    
                if num_specs>6
                    turb.ep(i,ind,iBeam) = (A(2)/0.53)^(3/2);
                    turb.N(i,ind,iBeam) = A(1);
                    turb.SNR(i,ind,iBeam) = A(2)./A(1);
                    turb.A(i,ind,iBeam) = A(2);
                    turb.corr(i,ind,iBeam) = corr_avg;
                    turb.spec(i,ind,:,iBeam) = permute(spec_fit,[2,3,1]);
                    turb.spec_num(i,ind,:,iBeam) = permute(avg_num,[2,3,1]);
                    turb.slope(i,ind,iBeam) = Coeffs_SL(2);
                    turb.N_slope(i,ind,iBeam) = Coeffs_SL(1);
                else
                    turb.ep(i,ind,iBeam) = NaN;
                    turb.N(i,ind,iBeam) = NaN;
                    turb.SNR(i,ind,iBeam) = NaN;
                    turb.A(i,ind,iBeam) = NaN;
                    turb.corr(i,ind,iBeam) = NaN;
                    turb.spec(i,ind,:,iBeam) = nan(1,1,length(ks));
                    turb.slope(i,ind,iBeam) = NaN;
                    turb.N_slope(i,ind,iBeam) = NaN;
                end
            end
        end
        turb.time(ind) = mean(time,'omitnan');
        ind = ind+1;
    end
end

turb.z = transpose(bin_centers);
turb.z_struct = transpose(z_grid_dec);  % depth grid for the structure-function estimates (ep_struct/A_struct/N_struct/struct_fun)
turb.k = transpose(ks);
turb.r = transpose(range_sub);
turb.beam_number = permute(HRbeams,[1,3,2]);

[turb.time,i] = unique(turb.time);
turb.ep = turb.ep(:,i,:);
turb.N = turb.N(:,i,:);
turb.SNR = turb.SNR(:,i,:);
turb.A = turb.A(:,i,:);
turb.ep_struct = turb.ep_struct(:,i,:);
turb.N_struct = turb.N_struct(:,i,:);
turb.A_struct = turb.A_struct(:,i,:);
turb.corr = turb.corr(:,i,:);
turb.spec = turb.spec(:,i,:,:);
turb.struct_fun = turb.struct_fun(:,i,:,:);
turb.slope = turb.slope(:,i,:);
turb.N_slope = turb.N_slope(:,i,:);
