clc
clear
close all

load('sine_xyz_fast.mat')
main_ekf_trilateration;
Xpo_trilat = Xpo;
clearvars -except Xpo_trilat
load('sine_xyz_fast.mat')
main_ekf;
Xpo_dist = Xpo;

%% Plot estimated altitude vs ground truth
figure
subplot(3,1,1)
plot(t,Xpo_dist(:,1),'r','Linewidth',2)
set(gca,'FontSize',16)
grid on
hold on
plot(t, Xpo_trilat(:,1),'g', 'LineWidth', 2)
plot(t_vicon,pos_vicon(:,1),'b','Linewidth',2)
xlabel('t [s]')
ylabel('x [m]')
legend('UWB Dist','UWB Trilat','VICON')

subplot(3,1,2)
plot(t,Xpo_dist(:,2),'r','Linewidth',2)
set(gca,'FontSize',16)
grid on
hold on
plot(t, Xpo_trilat(:,2),'g', 'LineWidth', 2)
plot(t_vicon,pos_vicon(:,2),'b','Linewidth',2)
xlabel('t [s]')
ylabel('y [m]')
legend('UWB Dist','UWB Trilat','VICON')

subplot(3,1,3)
plot(t,Xpo_dist(:,3),'r','Linewidth',2)
set(gca,'FontSize',16)
grid on
hold on
plot(t, Xpo_trilat(:,3),'g', 'LineWidth', 2)
plot(t_vicon,pos_vicon(:,3),'b','Linewidth',2)
xlabel('t [s]')
ylabel('z [m]')
legend('UWB Dist','UWB Trilat','VICON')
set(gcf,'color','w');

% Compute RMS errors in each direction
% Find closest ground truth data based on current time
for k = 1:K
    [~,idx_vicon(k)] = min(abs(t(k)-t_vicon));
end

% Compute error
rms_x_dist = rms(Xpo_dist(:,1) - pos_vicon(idx_vicon,1))
rms_y_dist = rms(Xpo_dist(:,2) - pos_vicon(idx_vicon,2))
rms_z_dist = rms(Xpo_dist(:,3) - pos_vicon(idx_vicon,3))

rms_x_trilat = rms(Xpo_trilat(:,1) - pos_vicon(idx_vicon,1))
rms_y_trilat = rms(Xpo_trilat(:,2) - pos_vicon(idx_vicon,2))
rms_z_trilat = rms(Xpo_trilat(:,3) - pos_vicon(idx_vicon,3))