function plot_result_adcp(WWmeta,Vel)
% code to take a quick look at the result
% need cbrewer colomap toolbox
% 
% Bofu Zheng

title_string = WWmeta.name_aqd;
title_string(strfind(title_string,'_'))=' ';

cb = cbrewer('div','RdBu',21);
figure
ax2 = subplot(211)
time = Vel{3};
dz = -Vel{4};
h = pcolor(time,dz,Vel{1});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.2 0.2])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
end
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' E-W velocity'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


ax1 = subplot(212)
h = pcolor(time,dz,Vel{2});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.2 0.2])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
end
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' N-S velocity'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)
linkaxes([ax1,ax2])

figure
ax1=subplot(211);
time = Vel{3};
dz = -Vel{4};
h = pcolor(time,dz,Vel{5});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.035 0.035])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
end
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' E-W shear'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


ax2=subplot(212);
h = pcolor(time,dz,Vel{6});
set(h, 'EdgeColor', 'none');
colormap(flipud(cb));
caxis([-0.035 0.035])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
enddatetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' N-S shear'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)
linkaxes([ax1,ax2])

figure
ax1=subplot(211);
time = Vel{3};
dz = -Vel{4};
h = pcolor(time,dz,nanmean(Vel{9},3));
set(h, 'EdgeColor', 'none');
colormap(ax1,'jet');
caxis([72 100])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
end
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' Relitive intensity'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)


ax2=subplot(212);
h = pcolor(time,dz,nanmean(Vel{10},3));
set(h, 'EdgeColor', 'none');
colormap(ax2,'jet');
caxis([0 10])
colorbar
danum = min(min(time));
if max(max(time))-min(min(time))>4
    set(gca,'xtick',danum:2:max(max(time)),'tickdir','in');
end
datetick('x','dd HH:MM','keepticks');
hylb = ylabel('Depth (m)');
set(hylb,'Color',[0.5 0.5 0.5],'FontSize',16)
xlim([min(min(time)) max(max(time))])
% ylim([0 25])
title(['\fontsize{15}',title_string,' intensity variance'])
hxlb = xlabel('Time (DD HH:MM)');
set(hxlb,'Color',[0.5 0.5 0.5],'FontSize',16)
linkaxes([ax1,ax2])

end
