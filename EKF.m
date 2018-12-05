function [Xpo,Ppo,t] = EKF(dataset,sensors)
%%%%%%%%%%%%%%%%%%%%%%%%%% DATASET FORMAT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% anchor_pos -> 3D location of anchors in VICON frame, ordered from 1-8
% flowdeck ->  first two columns correspond to the accumulated pixels in X, Y
%             third column corresponds to the Zranger measurement of height
% imu -> first three columns are accelerometer readings in xyz
%        last three columns are gyroscope measurements in xyz
% pos_vicon -> XYZ position measured by VICON system @ 200 Hz
% ref -> position setpoint sent to the drone
% uwb1 -> distance measurements to anchors 1-4
% uwb2 -> distance measurements to anchors 5-8

% INPUTS:
% dataset -> string containing the name of the dataset
% sensors -> struct with boolean variables defining which sensors
%            to use

% OUTPUTS:
% Xpo -> posterior estimate given by the EKF
% Ppo -> posterior covariance given by the EKF
% t -> timebase of the estimates, used for plotting

load(dataset)
%%%%%%%%%%%%%%%%%% INITIAL CONFIG %%%%%%%%%%%%%%%%%%%%%%%%%
% Decide which sensors we wish to fuse in the EKF
USE_IMU = sensors.use_imu;
USE_FLOW = sensors.use_flow;
USE_UWB = sensors.use_uwb;
USE_ZRANGER = sensors.use_zranger;
trilat = sensors.trilat;

% std deviations of initial states 
std_xy0 = 0.01; 
std_z0 = 0.01;
std_vel0 = 0.01;
std_rp0 = 0.01;
std_yaw0 = 0.01;

% Process noise 
w_accxy = 2;
w_accz = 1;
w_vel = 0;
w_pos = 0;
w_att = 0;
w_gyro_rp = 0.1;
w_gyro_yaw = 0.1;

% Constants
GRAVITY_MAGNITUDE = 9.81;
DEG_TO_RAD = pi/180.0;
e3 = [0; 0; 1];

% Flow deck constants
N_pix = 30.0; % number of pixels in the square image
theta_pix = 4.2 * DEG_TO_RAD; % angle of aperture
omega_factor = 1.25;

% Standard deviations of each sensor (tuning parameters)
std_flow = 0.1;
std_tof = 0.001; 
std_uwb = 0.03;
std_xy_trilat = 0.02;
std_z_trilat = 0.03;

%%%%%%%%%%%%%%%%% INITIALIZATION OF EKF %%%%%%%%%%%%%%%%%%%%

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

% Create a compound vector t with a sorted merge of all the time bases
time_bases = [];
if USE_IMU
    time_bases = [time_bases; t_imu];
end

if USE_UWB
    time_bases = [time_bases; t_uwb];
end

if USE_FLOW || USE_ZRANGER
    time_bases = [time_bases; t_flow];
end

t = unique(sort(time_bases));
K = length(t);

% Initial states/inputs
omega0 = zeros(3,1); % Initial angular velocity input

f_k = [0; 0; 0]; % Initial accelerometer input
f = zeros(K, 3);
f(1, 1:3) = f_k';

X0 = zeros(9,1); % Initial estimate for the state vector
X0(1) = 1;
R = eye(3); % Rotation from body frame to inertial frame

R_list = zeros(K, 3, 3);
R_list(1, :, :) = R;

% Initial posterior covariance
P0 = diag([std_xy0^2 std_xy0^2 std_z0^2,...
           std_vel0^2 std_vel0^2 std_vel0^2, ...
           std_rp0^2 std_rp0^2 std_yaw0^2]);

Xpr = zeros(K, 9);
Xpo = zeros(K, 9);
Xpr(1, :) = X0;
Xpo(1, :) = X0;

omega(1,:) = omega0;

Ppo = zeros(K, 9, 9);
Ppr = zeros(K, 9, 9);
Ppr(1,:,:) = P0;
Ppo(1,:,:) = P0;

%%%%%%%%%%%%%%%%%%% MAIN EKF LOOP %%%%%%%%%%%%%%%%%%%%%%%%%%
% Steps to follow:
% 1) At a given point in time, find out which measurements are available
% 2) Update the current states using said measurements, with the
%    corresponding EKF equations for the specific sensor
    
% Decide which anchors to fuse
anchor_id = [1,2,3,4,5,6,7,8];

for k = 2:K
        
    % Find what measurements are available at the current time
    imu_k = find(t_imu == t(k-1),1,'first');
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
    
    dt = t(k)-t(k-1);
    
    % Process noise
    Q = diag([(w_accxy*dt^2 + w_vel*dt + w_pos)^2,...
              (w_accxy*dt^2 + w_vel*dt + w_pos)^2,...
              (w_accz*dt^2 + w_vel*dt)^2,...
              (w_accxy*dt + w_vel)^2,...
              (w_accxy*dt + w_vel)^2,...
              (w_accz*dt + w_vel)^2,...
              (w_gyro_rp*dt + w_att)^2,...
              (w_gyro_rp*dt + w_att)^2,...
              (w_gyro_yaw*dt + w_att)^2]);

    if (~isempty(imu_k) && USE_IMU)
        % We have a new IMU measurement:
        % update the prior Xpr based on accelerometer and gyroscope data
        
        % Gyro measurements given in degrees, but EKF expects radians
        omega_k = imu(imu_k, 4:6)'*DEG_TO_RAD;        
        omega(k,:) = omega_k;
        
        Vpo = Xpo(k-1, 4:6)'; % Posterior estimate of velocity
        
        % Accelerometer measurements given in G's, but EKF expects m/s^2
        f_k = imu(imu_k, 1:3)'*GRAVITY_MAGNITUDE;
        f(k, 1:3) = f_k';
        
        d = omega_k*dt; % Attitude error
        
        % Motion Jacobian
        F = [ eye(3)    dt*R        -dt^2/2.0 * R * cross(f_k)
              zeros(3)  eye(3)      dt * R' * cross([0 0 GRAVITY_MAGNITUDE])
              zeros(3)  zeros(3)    expm( cross(d/2.0) )'];
        
        % Pass posterior covariance estimate through the motion model
        Ppr(k,:,:) = F * squeeze(Ppo(k-1,:,:)) * F' + Q;
        % Enforce symmetry
        Ppr(k,:,:) = 0.5*(squeeze(Ppr(k,:,:))+squeeze(Ppr(k,:,:))');
        
        % Orientation estimate
        R  = expm(cross(d))*R;
        Xpr(k, 7:9) = [0 0 0];
        
        R_list(k, :, :) = R;
        
%             dX = Vpo * dt + f_k*dt^2/2.0;
%             Xpr(k, 1:3) = Xpo(k-1, 1:3) + (R*dX)' - [0 0 GRAVITY_MAGNITUDE]*dt^2/2.0;
        Xpr(k, 1:3) = Xpo(k-1, 1:3) + Vpo'*dt  - [0 0 GRAVITY_MAGNITUDE]*dt^2/2.0;
        
        % Velocity prediction
        % Rotate gravity from along inertial z-axis into the body frame
        Xpr(k, 4:6) = Xpo(k-1, 4:6) + (f_k - GRAVITY_MAGNITUDE*R'*e3)'*dt;
        
        if (Xpr(k, 3) < 0)
            Xpr(k, 3:6) = [0 0 0 0];
        end     
               
    elseif(USE_IMU)
        % If we don't have IMU data, integrate the last acceleration
        % measurement (f_k) over one timestep to give our prior estimates
        Ppr(k,:,:) = squeeze(Ppo(k-1,:,:)) + Q;
        Ppr(k,:,:) = 0.5*(squeeze(Ppr(k-1,:,:))+squeeze(Ppr(k-1,:,:))');
        
        omega(k, :) = omega(k-1, :);
        
        % Position prediction
        dX = Xpo(k-1, 4:6)' * dt + f_k*dt^2/2.0;
        %Xpr(k, 1:3) = Xpo(k-1, 1:3) + (R*dX)' - [0 0 GRAVITY_MAGNITUDE]*dt^2/2.0;
        Xpr(k, 1:3) = Xpo(k-1, 1:3) + Xpo(k-1, 4:6)*dt  - [0 0 GRAVITY_MAGNITUDE]*dt^2/2.0;
        
        % Velocity prediction
        % Rotate gravity from along inertial z-axis to body frame
        Xpr(k, 4:6) = Xpo(k-1, 4:6) + (f_k - GRAVITY_MAGNITUDE*R'*e3)'*dt;
        
        f(k, 1:3) = f_k;
    else
        Ppr(k,:,:) = squeeze(Ppo(k-1,:,:)) + Q;
        Ppr(k,:,:) = 0.5*(squeeze(Ppr(k-1,:,:))+squeeze(Ppr(k-1,:,:))');
    end
    
    % Initially take our posterior estimates as the prior estimates
    % These are updated if we have flow deck or UWB measurements
    Xpo(k, :) = Xpr(k,:); 
    Ppo(k,:,:) = Ppr(k, :, :);
    updated = false;

    if(~isempty(flow_k) && (USE_FLOW || USE_ZRANGER))
                
        cos_alpha = R(3,3); % Cosine of the angle between body and inertial z-axes
        
        % Measurements from optical flow
        meas_nx = flowdeck(flow_k, 1);
        meas_ny = flowdeck(flow_k, 2);       

        % If measurements aren't outliers
        if abs(meas_nx) < 100 && abs(meas_ny) < 100
            
            v_x = Xpr(k, 4); % Velocity along body x-axis
            v_y = Xpr(k, 5); % Velocity along body y-axis

            % Saturate height to avoid singularities
            if (Xpr(k, 3) < 0.1)
                z = 0.1;
            else
                z = Xpr(k, 3);
            end

            omega_x = omega(k, 1);
            omega_y = omega(k, 2);

            pred_nx = dt*N_pix/theta_pix * (v_x  * cos_alpha/z - omega_factor*omega_y); 
            pred_ny = dt*N_pix/theta_pix * (v_y  * cos_alpha/z + omega_factor*omega_x);

            error_flow = [meas_nx - pred_nx;
                          meas_ny - pred_ny];

            % Covariance of flow deck measurements
            cov = diag([std_flow std_flow])^2; 

            % Jacobian for flow deck measurements
            h = zeros(2,9); 

            % Derivatives of pred_nx wrt z and v_x 
            h(1, 3) = -dt*N_pix/theta_pix * v_x  * cos_alpha / z^2;
            h(1, 4) = dt*N_pix/theta_pix * cos_alpha / z ;

            % Derivatives of pred_ny wrt z and v_y
            h(2, 3) = -dt*N_pix/theta_pix * v_y * cos_alpha / z^2;
            h(2, 5) = dt*N_pix/theta_pix * cos_alpha / z;

            [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),h,error_flow,cov);
            updated = true;

        end
                
        % Only fuse the z-ranger if cos(alpha) is large
        % If it's small, the body and inertial z-axes are nearly
        % perpendicular and thus the z-ranger will give poor measurements
        if (abs(cos_alpha) > 0.1 && cos_alpha > 0 && USE_ZRANGER)
            % Prediction and measurement for the z-ranger
            
            pred_dist = Xpo(k,3) / cos_alpha;
            meas_dist = flowdeck(flow_k,3);
            error_tof = meas_dist - pred_dist;
            var_tof = std_tof^2;
                       
            % Jacobian for z-ranger
            h = zeros(1,9); 
            
            % Derivative of pred_dist wrt z
            h(3) = 1 / cos_alpha;
            [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppo(k,:,:)),Xpo(k,:),h,error_tof,var_tof);

        end

    end
        
    if(~isempty(uwb_k) && USE_UWB)
        if (trilat)
            % We have a new UWB measurement from anchors 1-8
            % Do a trilateration to pinpoint 3D location of tag
            pos_trilat(:,uwb_k) = trilateration3D(anchor_pos(:,anchor_id),uwb(uwb_k,anchor_id)); 
            H = [eye(3,3) zeros(3,6)];
            var_xy = std_xy_trilat^2;
            var_z = std_z_trilat^2;
            Q = diag([var_xy var_xy var_z]);
            Err_uwb = pos_trilat(:,uwb_k) - Xpr(k,1:3)';
            [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),H,Err_uwb,Q);
        else
            % We have a new UWB measurement from anchors 1-8
            % update the states based on the distance measurement
            H = zeros(8,9);
            Q = std_uwb^2*eye(8,8);
            for i = 1:8
                dxi = Xpr(k,1) - anchor_pos(1,i);
                dyi = Xpr(k,2) - anchor_pos(2,i);
                dzi = Xpr(k,3) - anchor_pos(3,i);
                disti = sqrt(dxi^2 + dyi^2 + dzi^2);
                H(i,1:3) = [dxi/disti dyi/disti dzi/disti];
                Err_uwb(i) = uwb(uwb_k,i) - disti;
            end
            [Xpo(k,:),Ppo(k,:,:)] = update_state(squeeze(Ppr(k,:,:)),Xpr(k,:),H,Err_uwb',Q); 
        end
    end
    
    if (updated)
        v = Xpo(k, 7:9)*dt;     
        A = blkdiag(eye(3), eye(3), expm(cross(v/2.0)));
        Ppo(k,:,:) = A * squeeze(Ppo(k,:,:)) * A';
        Ppo(k,:,:) = 0.5*(squeeze(Ppo(k,:,:))+squeeze(Ppo(k,:,:))');
        R  = expm(cross(v))*R;
        Xpo(k, 7:9) = [0 0 0];
    end 
end

