clc
clear
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anchor_pos -> 3D location of anchors in VICON frame, ordered from 1-8
% flowdeck ->  first two columns correspond to the accumulated pixels in XY
%             third column corresponds to the Zranger measurement of height
% imu -> first three columns are accelerometer readings in xyz
%        last three columns are gyroscope measurements in xyz
% pos_vicon -> XYZ position measured by VICON system @ 200 Hz
% ref -> position setpoint sent to the drone
% uwb1 -> distance measurements to anchors 1-4
% uwb2 -> distance measurements to anchors 5-8


% Load the dataset
load('sine_xyz_fast_log30hz.mat')

% Plot followed trajectory in XYZ along with the commands sent
% FIGURE1: VICON ground truth and setpoints
figure(1)
subplot(3,1,1)
grid on;
hold on;
plot(t_vicon,pos_vicon(:,1),'b','Linewidth',1.5);
plot(t_cmds, ref(:,1),'--r', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('x [m]')
legend('VICON ground truth', 'Command')
set(gca,'FontSize',16)

subplot(3,1,2)
grid on;
hold on;
plot(t_vicon,pos_vicon(:,2),'b','Linewidth',1.5);
plot(t_cmds, ref(:,2),'--r', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('y [m]')
legend('VICON ground truth', 'Command')
set(gca,'FontSize',16)

subplot(3,1,3)
grid on;
hold on;
plot(t_vicon,pos_vicon(:,3),'b','Linewidth',1.5);
plot(t_cmds, ref(:,3),'--r', 'LineWidth', 1.5)
xlabel('Time [s]')
ylabel('z [m]')
legend('VICON ground truth', 'Command')
set(gca,'FontSize',16)