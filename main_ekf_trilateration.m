clear

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
load('sine_xyz_fast.mat')
%%%%%%%%%%%%%%%%%% INITIAL CONFIG %%%%%%%%%%%%%%%%%%%%%%%%%
% Decide which sensors we wish to fuse in the EKF
USE_IMU = false;
USE_FLOW = false;
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

% Create a unified time base for all UWB measurements
if length(t_uwb1) < length(t_uwb2)
    for k = 1:length(uwb1)
        [~,idx_uwb2(k)] = min(abs(t_uwb1(k)-t_uwb2));
    end
    uwb2 = uwb2(idx_uwb2,:);
    idx_uwb = find(abs(t_uwb1 - t_uwb2(idx_uwb2)) <= 0.033);
else
    for k = 1:length(uwb2)
        [~,idx_uwb1(k)] = min(abs(t_uwb2(k)-t_uwb1));
    end
    uwb1 = uwb1(idx_uwb1,:);
    idx_uwb = find(abs(t_uwb2 - t_uwb1(idx_uwb1)) <= 0.033);
end
uwb2 = uwb2(idx_uwb,:);
uwb1 = uwb1(idx_uwb,:);
t_uwb = t_uwb1(idx_uwb);
uwb = [uwb1 uwb2];
    
% Decide which anchors to fuse
anchor_id = [1,2,3,4,5,6,7,8];

% Create a compound vector t with a sorted merge of all the time bases
t = unique(sort([t_imu; t_uwb; t_flow]));
K = length(t);

% The EKF will provide an estimate of the state for each value of t
for k = 2:K
    % Find what measurements are available at the current time
    imu_k = find(t_imu == t(k-1));
    flow_k = find(t_flow == t(k-1),1,'first');
    uwb_k = find(t_uwb == t(k-1),1,'first');
    
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
    
    if(~isempty(uwb_k) && USE_UWB)
        % We have a new UWB measurement from anchors 1-8
        % Do a trilateration to pinpoint 3D location of tag
        pos(:,uwb_k) = trilateration3D(anchor_pos(:,anchor_id),uwb(uwb_k,anchor_id)); 
        H = [eye(3,3) zeros(3,6)];
        std_xy = 0.001;
        std_z = 0.008;
        Q = diag([std_xy std_xy std_z]);
        Err_uwb = pos(:,uwb_k) - Xpr(k,1:3)';
        [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),H,Err_uwb,Q);    
    end
end

%% Plot estimated altitude vs ground truth
figure
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
% 
% Compute RMS errors in each direction
% Find closest ground truth data based on current time
for k = 1:K
    [~,idx_vicon(k)] = min(abs(t(k)-t_vicon));
end
% 
% % Compute error
% rms_x = rms(Xpo(:,1) - pos_vicon(idx_vicon,1))
% rms_y = rms(Xpo(:,2) - pos_vicon(idx_vicon,2))
% rms_z = rms(Xpo(:,3) - pos_vicon(idx_vicon,3))

%% Visualize 3D
figure(2)
tk = 0:0.1:t(end);
pk = spline(t,pos_vicon(idx_vicon,:)',tk);
pk_est = spline(t,Xpo(:,1:3)',tk);
wnd = 20;
for k = 1:length(tk)
    [az,el] = view;
    for j = 1:8
        plot3(anchor_pos(1,j),anchor_pos(2,j),anchor_pos(3,j),'*k','markers',12)
        view([az,el]);
        hold on;
    end
    grid on
    if k <= wnd
        plot3(pk_est(1,1:k),pk_est(2,1:k),pk_est(3,1:k),'r','Linewidth',2)
        plot3(pk(1,1:k),pk(2,1:k),pk(3,1:k),'r','Linewidth',2)
    else
        plot3(pk_est(1,k-wnd:k),pk_est(2,k-wnd:k),pk_est(3,k-wnd:k),'r','Linewidth',2)
        plot3(pk(1,k-wnd:k),pk(2,k-wnd:k),pk(3,k-wnd:k),'b','Linewidth',2)
    end
    plot3(pk_est(1,k), pk_est(2,k),pk_est(3,k),...
                  'o','LineWidth',2,'MarkerEdgeColor','k',...
                      'MarkerFaceColor','r','markers',12);
    plot3(pk(1,k), pk(2,k),pk(3,k),...
                  'o','LineWidth',2,'MarkerEdgeColor','k',...
                      'MarkerFaceColor','b','markers',12); 
    xlabel('$x$ [m]','Interpreter','latex','Fontsize',16);
    ylabel('$y$ [m]','Interpreter','latex','Fontsize',16);
    zlabel('$z$ [m]','Interpreter','latex','Fontsize',16);
    set(gcf,'color','w');
    drawnow
    hold off
end

