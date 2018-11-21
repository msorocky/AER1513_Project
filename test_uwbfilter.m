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
outlier_th = 2;
% detecting outliers in UWB data
id = 4;
K = length(t_uwb2);
history_length = 4;
l =1;
for k=1:K
   if k > history_length
       mean_uwb = mean(uwb2(k-history_length:k-1,id));
       std_uwb = std(uwb2(k-history_length:k-1,id));
       dist = abs(mean_uwb - uwb2(k,id));
       if (dist < outlier_th*std_uwb)
           filtered(l) = uwb2(k,id);
           t_filtered(l) = t_uwb2(k);
           l = l + 1;
       end
   else
       filtered(l) =  uwb2(k,id);
       t_filtered(l) = t_uwb2(k);
       l = l + 1;
   end  
end
plot(t_filtered,filtered)
figure
plot(t_uwb2,uwb2(:,id))