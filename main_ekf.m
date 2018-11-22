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
% Select which mat file you want to load
% [FileName,PathName,FilterIndex] = uigetfile('.mat');
% load(FileName);
load('circle_xy.mat')

%%%%%%%%%%%%%%%%%% INITIAL CONFIG %%%%%%%%%%%%%%%%%%%%%%%%%
% Decide which sensors we wish to fuse in the EKF
USE_IMU = false;
USE_FLOW = true;
USE_UWB = true;

% std deviations of initial states 
std_xy0 = 1; 
std_z0 = 1;
std_vel0 = 0.01;
std_rp0 = 0.01;
std_yaw0 = 0.01;

% Process noise 
w_accxy = 2;
w_accz = 2;
w_vel = 0;
w_pos = 0;
w_att = 0;
w_gyro_rp = 0.1;
w_gyro_yaw = 0.1;

%%%%%%%%%%%%%%%%% INITIALIZATION OF EKF %%%%%%%%%%%%%%%%%%%%
X0 = zeros(9,1); % initial estimate for the state vector
X0(1) = 1;
R = eye(3); % attitude in rotation matrix form
q = [1;0;0;0]; % attitude in quaternion form 
P0 = diag([std_xy0^2 std_xy0^2 std_z0^2,...
           std_vel0^2 std_vel0^2 std_vel0^2, ...
           std_rp0^2 std_rp0^2 std_yaw0^2]);
Xpr(1,:) = X0;   % prior of estimated state
Xpo(1,:) = X0;   % posterior of estimated state
Ppo(1,:,:) = P0; % posterior covariance of estimated state
Ppr(1,:,:) = P0; % prior covariance of estimated state

%%%%%%%%%%%%%%%%%%% MAIN EKF LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps to follow:
% 1) At a given point in time, find out which measurements are available
% 2) Update the current states using said measurements, with the
%    corresponding EKF equations for the specific sensor

% Create a compound vector t with a sorted merge of all the time bases
t = unique(sort([t_imu; t_uwb1; t_uwb2; t_flow]));
K = length(t);

% The EKF will provide an estimate of the state for each value of t
for k = 2:K
    % Find what measurements are available at the current time
    imu_k = find(t_imu == t(k-1));
    flow_k = find(t_flow == t(k-1),1,'first');
    uwb1_k = find(t_uwb1 == t(k-1),1,'first');
    uwb2_k = find(t_uwb2 == t(k-1),1,'first');
    
    if(~isempty(imu_k) && USE_IMU)
        % We have a new IMU measurement:
        % update the prior Xpr based on accelerometer and gyroscope data
        
    else
        % No IMU data -> the prior knowledge of the state is just the
        % previous best estimate Xpo(k-1,:)
        Xpr(k,:) = Xpo(k-1,:);
    end
    
    % Compute Covariance prior Ppr
    dt = t(k)-t(k-1);
    Q = diag([(w_accxy*dt^2 + w_vel*dt + w_pos)^2,...
              (w_accxy*dt^2 + w_vel*dt + w_pos)^2,...
              (w_accz*dt^2 + w_vel*dt)^2,...
              (w_accxy*dt + w_vel)^2,...
              (w_accxy*dt + w_vel)^2,...
              (w_accz*dt + w_vel)^2,...
              (w_gyro_rp*dt + w_att)^2,...
              (w_gyro_rp*dt + w_att)^2,...
              (w_gyro_yaw*dt + w_att)^2]);
    Ppr(k,:,:) = squeeze(Ppo(k-1,:,:)) + Q;
    
    % this two lines are just for current testing, in every iteration of
    % the EKF both posteriors will always be updated, right now we force
    % this update to test single sensor updated, for example
    Ppo(k,:,:) = Ppr(k,:,:);
    Xpo(k,:) = Xpr(k,:); 

    if(~isempty(flow_k) && USE_FLOW)
        % We have a new Flowdeck measurement
        % update the states based in the measurement model
        pred_dist = Xpr(k,3) / R(3,3);
        meas_dist = flowdeck(flow_k,3);
        error = meas_dist - pred_dist;
        std_tof = 0.001;
        h = zeros(1,9);
        h(3) = 1 / R(3,3);
        [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),h,error,std_tof);    
    end
    
    if(~isempty(uwb1_k) && USE_UWB)
        % We have a new UWB measurement from anchors 1-4
        % update the states based in the measurement model
        H = zeros(4,9);
        std_uwb = 0.005;
        Q = diag([std_uwb std_uwb std_uwb std_uwb]);
        for i = 1:4
            dxi = Xpr(k,1) - anchor_pos(1,i);
            dyi = Xpr(k,2) - anchor_pos(2,i);
            dzi = Xpr(k,3) - anchor_pos(3,i);
            disti = sqrt(dxi^2 + dyi^2 + dzi^2);
            H(i,1:3) = [dxi/disti dyi/disti dzi/disti];
            Err_uwb1(i) = uwb1(uwb1_k,i) - disti;
        end
        [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),H,Err_uwb1',Q);    
    end
    
    if(~isempty(uwb2_k) && USE_UWB)
        % We have a new UWB measurement from anchors 5-8
        % update the states based in the measurement model
        H = zeros(4,9);
        std_uwb = 0.005;
        Q = diag([std_uwb std_uwb std_uwb std_uwb]);
        for i = 1:4
            dxi = Xpr(k,1) - anchor_pos(1,i+4);
            dyi = Xpr(k,2) - anchor_pos(2,i+4);
            dzi = Xpr(k,3) - anchor_pos(3,i+4);
            disti = sqrt(dxi^2 + dyi^2 + dzi^2);
            H(i,1:3) = [dxi/disti dyi/disti dzi/disti];
            Err_uwb2(i) = uwb2(uwb2_k,i) - disti;
        end
        [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),H,Err_uwb2',Q);
    end
end

%% Plot estimated altitude vs ground truth
figure(1)
subplot(3,1,1)
plot(t,Xpo(:,1),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,1),'b','Linewidth',2)
plot(t_cmds, ref(:,1),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('x [m]')
legend('Estimate','VICON ground truth', 'Command')

subplot(3,1,2)
plot(t,Xpo(:,2),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,2),'b','Linewidth',2)
plot(t_cmds, ref(:,2),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('y [m]')
legend('Estimate','VICON ground truth', 'Command')

subplot(3,1,3)
plot(t,Xpo(:,3),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,3),'b','Linewidth',2)
plot(t_cmds, ref(:,3),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('z [m]')
legend('Estimate','VICON ground truth', 'Command')