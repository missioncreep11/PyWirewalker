function out = WWcorr_beam(time,dpth,acce,pitch,roll,heading,vele_beam,saprate)
% function to calculate WW motion corrected velocity for down-ward looking
% ADCP on the WW
% saprate - sampling rate
% 
% out{1} - after motion correction (E-W)
% out{2} - after motion correction (N-S)
% out{3} - after motion correction (U-D)
% out{4} - ww motion (E-W)
% out{5} - ww motion (N-S)
% out{6} - ww motion (U-D)

% Bofu Zheng
% modified on Nov 2 2020


%% bandpass acceleration

% equal_depth_ENU_uncorr = Beam2ENU(squeeze(vele_beam(:,1,:)),squeeze(vele_beam(:,2,:)),...
%             squeeze(vele_beam(:,3,:)),squeeze(vele_beam(:,4,:)), pitch, roll, heading);
% vel_avg = nanmean(equal_depth_ENU_uncorr,3);

% dynamic accelerometer
accex = acce(:,1);   % get rid of tilt
accey = acce(:,2);
accez = acce(:,3);

% transfer acceleration from XYZ into ENU system
aENU = XYZ2ENU(accex, accey, accez, pitch, roll, heading);
aENU(:,3) = aENU(:,3)-1; % remove static acceleration

mask = abs(aENU)>0.3;mask = sum(mask,2)>0;
aENU(mask,:) = 0;

% filter acceleration
aENUf = [];  % creat a fake time series to get rid of the edge effect
aENUf(:,1) = [flip(aENU(:,1)); aENU(:,1); flip(aENU(:,1))];
aENUf(:,2) = [flip(aENU(:,2)); aENU(:,2); flip(aENU(:,2))];
aENUf(:,3) = [flip(aENU(:,3)); aENU(:,3); flip(aENU(:,3))];

% [b,a] = butter(1,[0.1 1.2]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
[b,a] = butter(1,[0.1 1.2]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
acu = filtfilt(b,a,double(aENUf(:,1)));
acu = detrend(acu,'constant');  % get rid of undesired trend
acu = acu(length(acu)*1/3+1:length(acu)*2/3);  % chunk back to original length
acu = acu - nanmean(acu); % demean

acv = filtfilt(b,a,double(aENUf(:,2)));
acv = detrend(acv,'constant');
acv = acv(length(acv)*1/3+1:length(acv)*2/3);
acv = acv - nanmean(acv); % demean

[b,a] = butter(1,[0.1 1.2]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
% [b,a] = butter(1,[0.01 3]/(saprate/2),'bandpass');  % set up filter, 0.1-1.2 Hz
acw = filtfilt(b,a,double(aENUf(:,3)));
acw = detrend(acw,'constant');
acw = acw(length(acw)*1/3+1:length(acw)*2/3);
acw = acw - nanmean(acw); % demean


dpw = (diff(dpth(1:end-1))-flipud(diff(flipud(dpth(2:end)))))/(2/saprate);
dpw = [flip(dpw); dpw; flip(dpw)];
[b,a] = butter(1,0.3/(saprate/2),'low');  % set up filter, lowpass cutoff 0.1Hz
dpw = filtfilt(b,a,double(dpw));
dpw = dpw(length(dpw)*1/3:length(dpw)*2/3+1);


%% calculate WW motion - translational velocity
WWu = zeros(size(time));
WWv = zeros(size(time));
WWw = zeros(size(time));

for i = 1:length(WWu)-1
    WWu(i+1) = WWu(i) + acu(i)*9.81*86400*(time(i+1)-time(i));  % acceleration in unit of g
    WWv(i+1) = WWv(i) + acv(i)*9.81*86400*(time(i+1)-time(i));
    WWw(i+1) = WWw(i) + acw(i)*9.81*86400*(time(i+1)-time(i));
end

% [f,coh_pres,phase_pres,~,xx_list_pres,yy_list_pres,xy_list_pres]=cohWindowed(vel_avg(2:651,3),dpw(1:650),76,0.5,'hamming');
% [f,coh_imu,phase_imu,~,xx_list_imu,yy_list_imu,xy_list_imu]=cohWindowed(vel_avg(2:651,3),-WWw(2:651),76,0.5,'hamming');
% [f,coh_comb,phase_comb,~,xx_list_comb,yy_list_comb,xy_list_comb]=cohWindowed(vel_avg(2:651,3),dpw(1:650)-WWw(2:651),76,0.5,'hamming');
% 
% Gain_pres = sqrt(xy_list_pres(1:length(f)).^2./xx_list_pres(1:length(f)).^2);
% Gain_imu = sqrt(xy_list_imu(1:length(f)).^2./xx_list_imu(1:length(f)).^2);
% Gain_comb = sqrt(xy_list_comb(1:length(f)).^2./xx_list_comb(1:length(f)).^2);
% 
% figure;subplot(131);semilogx(f*8,coh_pres(1:length(f)))
% hold on;semilogx(f*8,coh_imu(1:length(f)))
% hold on;semilogx(f*8,coh_comb(1:length(f)),'k','LineWidth',2)
% 
% subplot(132);semilogx(f*8,Gain_pres(1:length(f)))
% hold on;semilogx(f*8,Gain_imu(1:length(f)))
% hold on;semilogx(f*8,Gain_comb(1:length(f)),'k','LineWidth',2)
% 
% subplot(133);semilogx(f*8,phase_pres(1:length(f)))
% hold on;semilogx(f*8,phase_imu(1:length(f)))
% hold on;semilogx(f*8,phase_comb(1:length(f)),'k','LineWidth',2)

%% motion corrected velocity

vel_xyz = XYZ2ENU(WWu,WWv,-dpw+WWw,pitch,roll,heading,'reverse');
% vel_xyz = XYZ2ENU(WWu,WWv,WWw,pitch,roll,heading,'reverse');
% Signuature 1000 beam geometry
phi = [65,65,65,65]*pi/180;
azi = [0,-90,180,90]*pi/180;

% Get signature 1000 beam unit vectors
clear bX bY bZ
for ibeam = 1:4
    [bX(ibeam),bY(ibeam),bZ(ibeam)] = GetUnitVectors(phi(ibeam),azi(ibeam),0,0);
end
        
vel_corr_beam = zeros(size(vele_beam,1),size(vele_beam,2));
for ibeam = 1:4
    vel_corr_beam(:,ibeam) = vel_xyz*[bX(ibeam);bY(ibeam);bZ(ibeam)];
end

vel_beam_mc = vele_beam+vel_corr_beam;

%%
out{1} = vel_beam_mc;  % motion corrected
out{2} = vel_corr_beam;
out{3} = WWu;    % WW motion
out{4} = WWv;
out{5} = -dpw+WWw;
% out{5} = WWw;

end