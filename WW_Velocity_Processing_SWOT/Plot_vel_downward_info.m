depth_cut = zeros(length(time),4);
for b = 1:4
    for n = 1:length(time)
        depths = dpth_temp(n,:);
        depth_cut(n,b) = depths(cut_ind(n,b));
    end
end

figure;h=pcolor(time_temp,dpth_temp,velemc{1});h.LineStyle='none';

figure;h=pcolor(time_temp,dpth_temp,squeeze(beam_vel(:,1,:)));h.LineStyle='none';

figure;
subplot(2,2,1)
h=pcolor(time_temp,dpth_temp,squeeze(amp(:,1,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,1),'k','LineWidth',2)
subplot(2,2,2)
h=pcolor(time_temp,dpth_temp,squeeze(amp(:,2,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,2),'k','LineWidth',2)
subplot(2,2,3)
h=pcolor(time_temp,dpth_temp,squeeze(amp(:,3,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,3),'k','LineWidth',2)
subplot(2,2,4)
h=pcolor(time_temp,dpth_temp,squeeze(amp(:,4,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,4),'k','LineWidth',2)

figure;
subplot(2,2,1)
h=pcolor(time_temp(:,2:end),dpth_temp(:,2:end),squeeze(diffamp(:,1,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,1),'k','LineWidth',2)
subplot(2,2,2)
h=pcolor(time_temp(:,2:end),dpth_temp(:,2:end),squeeze(diffamp(:,2,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,2),'k','LineWidth',2)
subplot(2,2,3)
h=pcolor(time_temp(:,2:end),dpth_temp(:,2:end),squeeze(diffamp(:,3,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,3),'k','LineWidth',2)
subplot(2,2,4)
h=pcolor(time_temp(:,2:end),dpth_temp(:,2:end),squeeze(diffamp(:,4,:)));h.LineStyle='none';
hold on;plot(time_temp(:,1),depth_cut(:,4),'k','LineWidth',2)
