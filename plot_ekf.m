clc
close all
clear

% Plot EKFs performance given a dataset and the use of certain sensors
dataset = 'sine_xyz_fast.mat';

% Create all the trials you want to test

% Only trilateration
sensors.use_imu = false;
sensors.use_flow = false;
sensors.use_uwb = true;
sensors.use_zranger = false;
sensors.trilat = true;
sensors_choice(1) = sensors;
name{1} = 'Trilateration only';

% IMU + FLOW + ZRANGER 
sensors.use_imu = true;
sensors.use_flow = true;
sensors.use_uwb = false;
sensors.use_zranger = true;
sensors.trilat = false;
sensors_choice(2) = sensors;
name{2} = 'IMU/Flowdeck';

% All Sensors
sensors.use_imu = true;
sensors.use_flow = true;
sensors.use_uwb = true;
sensors.use_zranger = true;
sensors.trilat = false;
sensors_choice(3) = sensors;
name{3} = 'All Sensors';

% Load reference signal and groundtruth from the dataset
load(dataset,'anchor_pos','ref','t_cmds','pos_vicon','t_vicon');

for i = 1:length(sensors_choice)
    % Run EKF
    [Xpo{i},Ppo{i},t{i}] = EKF(dataset,sensors_choice(i));
    
    % Find closest ground truth data based on current time
    for k = 1:length(t{i})
        [~,idx_vicon{i}(k)] = min(abs(t{i}(k)-t_vicon));
    end
    
    % Compute RMS error between estimate and groundtruth data
    rms_x = rms(Xpo{i}(:,1) - pos_vicon(idx_vicon{i},1));
    rms_y = rms(Xpo{i}(:,2) - pos_vicon(idx_vicon{i},2));
    rms_z = rms(Xpo{i}(:,3) - pos_vicon(idx_vicon{i},3));

    fprintf("RMS error in X direction for " + name{i} + " --> %.2f cm\n",rms_x*100);
    fprintf("RMS error in Y direction for " + name{i} + " --> %.2f cm\n",rms_y*100);
    fprintf("RMS error in Z direction for " + name{i} + " --> %.2f cm\n",rms_z*100);
    fprintf("------------------------------------------\n")
end

%% Plot estimated position vs ground truth
figure(1)
colors = get(gca,'colororder');

for i = 1:length(sensors_choice)
    subplot(3,1,1)
    grid on
    hold on
    h_plot(i) = plot(t{i},Xpo{i}(:,1),'Color',colors(i+1,:),'Linewidth',2);
    set(gca,'FontSize',16)
    h_label{i} = [name{i}];
    h_vic = plot(t_vicon,pos_vicon(:,1),'b','Linewidth',2);
    xlabel('t [s]')
    ylabel('x [m]')
    if i == length(sensors_choice)
        legend([h_plot h_vic], [h_label, 'VICON'] )
    end

    subplot(3,1,2)
    h_plot2(i) = plot(t{i},Xpo{i}(:,2),'Color',colors(i+1,:),'Linewidth',2);
    set(gca,'FontSize',16)
    h_label2{i} = [name{i}];
    grid on
    hold on
    plot(t_vicon,pos_vicon(:,2),'b','Linewidth',2)
    xlabel('t [s]')
    ylabel('y [m]')
    if i == length(sensors_choice)
        legend([h_plot2 h_vic], [h_label2, 'VICON'] )
    end
    
    subplot(3,1,3)
    h_plot3(i) = plot(t{i},Xpo{i}(:,3),'Color',colors(i+1,:),'Linewidth',2);
    set(gca,'FontSize',16)
    h_label3{i} = [name{i}];
    grid on
    hold on
    plot(t_vicon,pos_vicon(:,3),'b','Linewidth',2)
    xlabel('t [s]')
    ylabel('z [m]')
    if i == length(sensors_choice)
        legend([h_plot3 h_vic], [h_label3, 'VICON'] )
    end
end
set(gcf,'color','w');

%% Visualize 3D
figure(2)
set(gcf, 'Position', get(0, 'Screensize'));
view([30,18]);
% Create a unified time base for all estimates
for i = 1:length(sensors_choice)
    length_t(i) = length(t{i});
end

[min_length,idx] = min(length_t);

for i = 1:length(sensors_choice)
   if i ~= idx
        % Find closest ground truth data based on current time
        for k = 1:length(t{idx})
            [~,idx_i{i}(k)] = min(abs(t{idx}(k)-t{i}));
        end
   end
end

tk = 0:0.1:t{idx}(end);
pk_vicon = spline(t{idx},pos_vicon(idx_vicon{idx},:)',tk);
for i = 1:length(sensors_choice)
    if i ~= idx
        pk_est{i} = spline(t{idx},Xpo{i}(idx_i{i},1:3)',tk);
    else
        pk_est{idx} = spline(t{idx},Xpo{idx}(:,1:3)',tk);
    end
end

wnd = 20;
for k = 1:length(tk)
    [az,el] = view;
    for j = 1:8
        h_anchor = plot3(anchor_pos(1,j),anchor_pos(2,j),anchor_pos(3,j),'*r','markers',12);
        set(gca,'FontSize',16)
        view([az,el]);
        hold on;
    end
    grid on
    if k <= wnd
        for i = 1:length(sensors_choice) 
            h_plot(i) = plot3(pk_est{i}(1,1:k),pk_est{i}(2,1:k),pk_est{i}(3,1:k),'Color',colors(i+1,:),'Linewidth',2);
            h_label{i} = [name{i}];
            h_vic = plot3(pk_vicon(1,1:k),pk_vicon(2,1:k),pk_vicon(3,1:k),'b','Linewidth',2);
        end
    else
        for i = 1:length(sensors_choice) 
            h_plot(i) = plot3(pk_est{i}(1,k-wnd:k),pk_est{i}(2,k-wnd:k),pk_est{i}(3,k-wnd:k),'Color',colors(i+1,:),'Linewidth',2);
            h_label{i} = [name{i}];
            h_vic = plot3(pk_vicon(1,k-wnd:k),pk_vicon(2,k-wnd:k),pk_vicon(3,k-wnd:k),'b','Linewidth',2);
        end
    end
    for i = 1:length(sensors_choice)
        plot3(pk_est{i}(1,k), pk_est{i}(2,k),pk_est{i}(3,k),...
                      'o','LineWidth',2,'MarkerEdgeColor','k',...
                          'MarkerFaceColor',colors(i+1,:),'markers',12);
        plot3(pk_vicon(1,k), pk_vicon(2,k), pk_vicon(3,k),...
                      'o','LineWidth',2,'MarkerEdgeColor','k',...
                          'MarkerFaceColor','b','markers',12); 
    end
                  
     
    xlabel('x [m]');
    ylabel('y [m]');
    zlabel('z [m]');
    xlim([-6 6])
    ylim([-6 6])
    zlim([0 3.5])
    legend([h_anchor h_plot h_vic], ['UWB Nodes', h_label, 'VICON'],'AutoUpdate','off')
    set(gcf,'color','w');
    drawnow
    hold off
end