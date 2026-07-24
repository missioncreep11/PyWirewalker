% plot_uv_tke_profile.m
% Plot a single Wirewalker ADCP profile:
%   left panel  : u (East) and v (North) velocity vs depth
%   right panel : TKE dissipation rate epsilon vs depth, on a log x-axis
%
% You enter the profile number at the prompt. The velocity product (ADCP) and
% the turbulence product (turb) are gridded independently, so the turbulence
% profile is matched to the chosen velocity cast by NEAREST TIME.

clear; clc;

%% ---- settings (match process_WW_ADCP_main_NOPP_d2_11.m) ----
root       = '/Users/drew/PyWirewalker/WW_Velocity_Processing_SWOT';
Wirewalker = 'WW';
Deployment = 'd2';
name_aqd   = 'S100430A038_NOPP_d2';
splitnum   = 1;          % which Grid chunk to read

gridpath = fullfile(root,'proc','Grid',Wirewalker,Deployment);
velfile  = fullfile(gridpath, sprintf('%s_%d.mat', name_aqd, splitnum));
turbfile = fullfile(gridpath, sprintf('%s_%d_HR_Turbulence.mat', name_aqd, splitnum));

%% ---- load ----
S = load(velfile,'ADCP'); ADCP = S.ADCP;
nprof = size(ADCP.velE,2);

hasturb = isfile(turbfile);
if hasturb
    Tt = load(turbfile,'turb'); turb = Tt.turb;
end

%% ---- choose profile ----
pnum = input(sprintf('Enter profile number (1-%d): ', nprof));
if isempty(pnum) || pnum < 1 || pnum > nprof || pnum ~= round(pnum)
    error('Profile number must be an integer between 1 and %d.', nprof);
end

u    = ADCP.velE(:,pnum);
v    = ADCP.velN(:,pnum);
zvel = ADCP.dz(:,pnum);                       % depth (m), positive downward
tvel = mean(ADCP.time(:,pnum),'omitnan');     % cast time (datenum)
zbot = max(zvel(~isnan(u) | ~isnan(v)));
if isempty(zbot) || ~isfinite(zbot); zbot = max(zvel); end

%% ---- figure ----
figure('Color','w','Position',[100 100 850 700]);

% --- left: u, v ---
ax1 = subplot(1,2,1);
plot(u, zvel, '-', 'Color',[0 0.45 0.74], 'LineWidth',1.5); hold on
plot(v, zvel, '-', 'Color',[0.85 0.33 0.10],'LineWidth',1.5);
xline(0,'k:');
set(ax1,'YDir','reverse'); grid on; box on
xlabel('velocity (m s^{-1})'); ylabel('depth (m)');
legend({'u (East)','v (North)'},'Location','best');
title(sprintf('Profile %d   %s', pnum, datestr(tvel,'yyyy-mm-dd HH:MM:SS')));

% --- right: epsilon (log x) ---
ax2 = subplot(1,2,2);
if hasturb && isfield(turb,'ep') && ~isempty(turb.ep)
    [dt, it] = min(abs(turb.time(:) - tvel));   % nearest turbulence profile in time
    ep = squeeze(turb.ep(:,it,:));              % (depth x beam)
    semilogx(ep, turb.z, '-', 'LineWidth',1.5);
    set(ax2,'YDir','reverse'); grid on; box on
    xlabel('\epsilon (W kg^{-1})');
    title(sprintf('TKE dissipation (turb prof %d, \\Deltat = %.0f s)', it, dt*86400));
    if size(ep,2) > 1
        legend(arrayfun(@(b) sprintf('beam %d',b), turb.beam_number(:), 'uni',0), ...
               'Location','best');
    end
else
    text(0.5,0.5,'no turbulence product found', 'Units','normalized', ...
         'HorizontalAlignment','center');
    set(ax2,'YDir','reverse'); box on
    xlabel('\epsilon (W kg^{-1})');
end

linkaxes([ax1 ax2],'y');
ylim(ax1,[0 zbot]);

