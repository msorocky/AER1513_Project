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

% Decide which sensors we wish to fuse in the EKF
VISUALIZE_UWB1 = true;
VISUALIZE_UWB2 = false;

% Create a compound vector t with a sorted merge of all the time bases
t = unique(sort([t_imu; t_uwb1; t_uwb2; t_flow]));
K = length(t);

figure(1)
% The EKF will provide an estimate of the state for each value of t
for k = 2:K
    % Find what measurements are available at the current time
    uwb1_k = find(t_uwb1 == t(k-1),1,'first');
    uwb2_k = find(t_uwb2 == t(k-1),1,'first');
    
    % Find closest ground truth data based on current time
    [~,idx_vicon] = min(abs(t(k)-t_vicon));
    [az,el] = view;
    if(~isempty(uwb1_k) && VISUALIZE_UWB1)
        % We have a new UWB measurement from anchors 1-4
        % update the states based in the measurement model
        CMap = parula(4);
        A = [[anchor_pos(:,1)' uwb1(uwb1_k,1)];
             [anchor_pos(:,2)' uwb1(uwb1_k,2)];
             [anchor_pos(:,3)' uwb1(uwb1_k,3)];
             [anchor_pos(:,4)' uwb1(uwb1_k,4)]];
        [Xi, Yi, Zi] = sphere;
        
        for j = 1:4
          XX = Xi * A(j, 4) + A(j, 1);
          YY = Yi * A(j, 4) + A(j, 2);
          ZZ = Zi * A(j, 4) + A(j, 3);
          plot3(anchor_pos(1,j),anchor_pos(2,j),anchor_pos(3,j),'.b','markers',15)
          view([az,el]);
          hold on;
          surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(j, :), ...
                  'EdgeColor', [0.4, 0.4, 0.4]); 
        end
        plot3(pos_vicon(idx_vicon,1),pos_vicon(idx_vicon,2),pos_vicon(idx_vicon,3),'.r','markers',50)
        grid on;
        xlim([-10,10])
        ylim([-10,10])
        zlim([0,10]);
        xlabel('$x$ [m]','Interpreter','latex','Fontsize',16);
        ylabel('$y$ [m]','Interpreter','latex','Fontsize',16);
        zlabel('$z$ [m]','Interpreter','latex','Fontsize',16);
        set(gcf,'color','w');
        drawnow
        hold off
    end
    
    if(~isempty(uwb2_k) && VISUALIZE_UWB2)
        % We have a new UWB measurement from anchors 5-8
        % update the states based in the measurement model
        CMap = parula(4);
        A = [[anchor_pos(:,5)' uwb2(uwb2_k,1)];
             [anchor_pos(:,6)' uwb2(uwb2_k,2)];
             [anchor_pos(:,7)' uwb2(uwb2_k,3)];
             [anchor_pos(:,8)' uwb2(uwb2_k,4)]];
        [Xi, Yi, Zi] = sphere;
        
        for j = 1:4
          XX = Xi * A(j, 4) + A(j, 1);
          YY = Yi * A(j, 4) + A(j, 2);
          ZZ = Zi * A(j, 4) + A(j, 3);
          plot3(anchor_pos(1,j+4),anchor_pos(2,j+4),anchor_pos(3,j+4),'.b','markers',15)
          view([az,el]);
          hold on;
          surface(XX, YY, ZZ, 'FaceAlpha', 0.5, 'FaceColor', CMap(j, :), ...
                  'EdgeColor', [0.4, 0.4, 0.4]); 
        end
        plot3(pos_vicon(idx_vicon,1),pos_vicon(idx_vicon,2),pos_vicon(idx_vicon,3),'.r','markers',50)
        grid on;
        xlim([-10,10])
        ylim([-10,10])
        zlim([0,10]);
        xlabel('$x$ [m]','Interpreter','latex','Fontsize',16);
        ylabel('$y$ [m]','Interpreter','latex','Fontsize',16);
        zlabel('$z$ [m]','Interpreter','latex','Fontsize',16);
        drawnow
        hold off
        
    end
end
