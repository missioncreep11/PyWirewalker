% main code to process WW ADCP data
% ADCP is in downward looking configuration
%
% Bofu Zheng/ Drew Lucas/ Arnaud Le Boyer/ 
% Oct.11 2021
% boz080@ucsd.edu

% Modified for use by Caeli Griffin April 2025

%% set path
% here the code is fed with .mat files of the measured velocity. To obtain
% .mat files, user needs to run Signature Deployment software first to
% convert .ad2cp files into .mat files

clc; clear all;
MainPath = '/Users/caeligriffin/Desktop/swot_nopp_s2/code/START202305_END202306_DEPLOY1-SUBSET/proc/';
% MainPath = '/Users/Devon/Documents/GradSchool/Nortek_Turbulence/TLC_25/';
Wirewalker = 'WW';   
Deployment = 'D1-SUB';

% Path to raw data
WWmeta.aqdpath=['/Users/caeligriffin/Desktop/swot_nopp_s2/code/START202305_END202306_DEPLOY1-SUBSET/raw/'];
% WWmeta.aqdpath=['/Users/Devon/Documents/GradSchool/Nortek_Turbulence/TLC_25/raw/'];
% root for WW_ADCP toolbox
WWmeta.root_script='/Users/caeligriffin/Desktop/swot_nopp_s2/WW_Velocity_Processing_SWOT_CLEAN'; 
% WWmeta.root_script='/Users/Devon/Documents/GradSchool/wirewalker/WW_Velocity_Processing_SWOT/';
% Name of the processed data, can be changed according to different cruises
WWmeta.name_aqd=['SN23_' Wirewalker '_' Deployment]; 

WWmeta = SetupPath(WWmeta,MainPath,Wirewalker,Deployment);

WWmeta % display what has been entered
cd(WWmeta.root_script) % change directory to the location...
dd0 = dir([WWmeta.aqdpath '*.mat']);
WWmeta.dd0 = dd0;
dd0 = dd0(~startsWith({dd0.name}, '._')); % Added by Caeli Griffin to remove files that start with '._'.

%% set variables
% adjustable variables include:
% 
% NUM_combining_files: number of .mat files to be combined as a group from 
%                      raw output of the ADCP, typically is set to be 20.
%                      if this number is too large, combined file is too 
%                      big to be saved
% blockdis           : blocking distance, can be found in the Config
%                      structure of the raw .mat file
% cellsize           : cell size, can be found in the Config structure of
%                      the raw .mat file
% saprate            : sampling rate, in Hz, can be found in the Config
%                      structure of the raw .mat file
% boxsize            : vertical range for averaging 
%                      - determining the vertical resolution of the final product
%                      typically is set to be the same as cell size
% z_max              : max depth of the WW profile, positive value
% k                  : determine whether to process downcast data. if k==1,
%                      downcast data will be saved. if k~=1, only upcast
%                      data will be processed and saved.
% thhold             : threshold value to determine whether it is too short
%                      to be a profile
variables.NUM_combining_files = 1;  % need to specify the number for combining files, the default is: 1
variables.blockdis = 0.1;            % need to specify the blanking distance, the default is: 0.1m
variables.cellsize = 0.5;             % need to specify the cell size, the default is: 1m
variables.saprate = 8;              % need to specify the sampling rate, the default is: 16Hz
variables.boxsize = 1;             % typically set to </= 1m (instructions updated by Caeli Griffin, following conversations with Drew Lucas & Devon Northcott)
variables.z_max   = 510;              % 500m profile, the default is: 500m, should slightly exceed actually max depth 
variables.k = 0;                     % 1 or not 1 (Is upcast data processed)
variables.thhold = 2;              % need to specify the number (2)
variables.direction = 'up';        % up or down, up or down facing adcp?
variables.sail_corr = 1;           % Correct for horizontal motion of the wirewalker? yes = 1, no = 0
variables.z_unit = [0,0,1]; % [-1/sqrt(2)*sind(22.5),-1/sqrt(2)*sind(22.5),cosd(22.5)]; % Unit vector of Nortek z-axis relative to wirewalker z-axis

variables.HRturb = 0;            % Process HR mode data for turbulence? yes = 1, no = 0

if variables.HRturb == 1
    variables.HRbeams = [1,2,3,4];  % Beams with HR mode enabled
    variables.HRblockdis = 0.1;     % blocking distance with HR mode
    variables.HRcellsize = 0.06;    % cell size for HR measurements 
    variables.HRboxsize = 50*variables.HRcellsize; % Depth resolution of final turbulence estimates
end

%% sort files (slow!)
% this is to make sure raw .mat files are in the right order
WWmeta = sort_file(WWmeta)

%% combine separate raw .mat files together and then chunk into profiles
for q = 1:variables.NUM_combining_files:length(dd0)  
    if q+variables.NUM_combining_files-1>length(dd0)
        num = length(dd0)-q+1;
    else
        num = variables.NUM_combining_files;
    end
    
    merge_signature(WWmeta,q,num);    % merge separate .mat files from ADCP output and then form a group
    create_profiles(WWmeta,q,num,variables.thhold,variables.k);  % will chunk ww profiles into upcast and downcast, and save them separately
    disp(['current file location: ',num2str(q),'_',num2str(q+num-1)])  % show where we are
end
disp('identify profiles: finished')

%% combine cut-off profiles
% there may be some profiles (specifically the first profile or last profile in a group)
% with first half in the previous group and second half in the current group. 
% Therefore, we are going to combine the cut-off profiles. 
copyfile([WWmeta.propath,'*.mat'],WWmeta.propath_rearrange)  % copy file from the old folder to the new folder
disp('copying file: finished')
combine_cutoff(WWmeta,variables.NUM_combining_files,variables.k)  % combine_cutoff is performed in the new folder
disp('combining: finished')

%% WWvel analysis
% here the motion correction and box averaging are performed

WWmeta.dd0 = dd0;
dd0 = dd0(~startsWith({dd0.name}, '._')); % Added by Caeli Griffin to remove files that start with '._'.

splitfiles = 25;splitnum=1;
numfiles = ceil(length(WWmeta.dd0)/variables.NUM_combining_files);

while splitfiles*(splitnum-1)<numfiles
    Vel = WWvel_upward(WWmeta,variables,0,splitfiles,splitnum);  % function to generate estimated velocity field

    ADCP.time = Vel{3};
    ADCP.dz   = Vel{4};
    ADCP.velE = Vel{1};
    ADCP.velN = Vel{2};
    ADCP.velU = Vel{14};
    ADCP.shearE = Vel{5};
    ADCP.shearN = Vel{6};
    ADCP.surf_vel = Vel{7};
    ADCP.Nav = Vel{8};
    ADCP.amp = Vel{9};
    ADCP.amp_var = Vel{10};
    ADCP.velE_var = Vel{12};
    ADCP.velN_var = Vel{13};
    ADCP.velU_var = Vel{15};
    ADCP.N = Vel{11};
    if variables.sail_corr == 1
        ADCP.velE_corr = Vel{16};
        ADCP.velN_corr = Vel{17};
    end

    ADCP.Notes = '';

    save([WWmeta.gridpath,WWmeta.name_aqd,'_' num2str(splitnum) '.mat'],'ADCP');  % save result

    if variables.HRturb==1
        turb = WWturb_upward(WWmeta,variables,splitfiles,splitnum);
        save([WWmeta.gridpath,WWmeta.name_aqd,'_' num2str(splitnum) '_HR_Turbulence.mat'],'turb');
    end
    
    splitnum = splitnum+1;
end

%% take a quick look at the result
% plot_result_adcp(WWmeta,Vel)




