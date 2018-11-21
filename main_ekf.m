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
load('sine_xyz.mat')

%%%%%%%%%%%%%%%%%% INITIAL CONFIG %%%%%%%%%%%%%%%%%%%%%%%%%
% Decide which sensors we wish to fuse in the EKF
USE_IMU = true;
USE_FLOW = true;
USE_UWB = true;

% MAIN EKF LOOP
% Steps to follow:
% 1) At a given point in time, find out which measurements are available
% 2) Update the current states using said measurements, with the
%    corresponding EKF equations for the specific sensor

% Create a compound vector t with a sorted merge of all the time bases
t = unique(sort([t_imu; t_uwb1; t_uwb2; t_flow]));
K = length(t);

% The EKF will provide an estimate of the state for each value of t
for k = 1:K
    % Find what measurements are available at the current time
    imu_k = find(t_imu == t(k));
    flow_k = find(t_flow == t(k));
    uwb1_k = find(t_uwb1 == t(k));
    uwb2_k = find(t_uwb2 == t(k));
    
    if(~isempty(imu_k) && USE_IMU)
        % We have a new IMU measurement:
        % update your state based on the prediction model
    end
    
    if(~isempty(flow_k) && USE_FLOW)
        % We have a new Flowdeck measurement
        % update the states based in the measurement model
    end
    
    if(~isempty(uwb1_k) && USE_UWB)
        % We have a new UWB measurement from anchors 1-4
        % update the states based in the measurement model
    end
    
    if(~isempty(uwb2_k) && USE_UWB)
        % We have a new UWB measurement from anchors 5-8
        % update the states based in the measurement model
    end
end





