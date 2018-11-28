clc
close all
clear 

load('sine_xyz_log30hz.mat')

% Find closest ground truth data based on current time
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

for k = 1:length(t_uwb)
   pos(:,k) = trilateration3D(anchor_pos(:,anchor_id),uwb(k,anchor_id)); 
end

%% Plot estimated altitude vs ground truth
figure(1)
subplot(3,1,1)
plot(t_uwb,pos(1,:),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,1),'b','Linewidth',2)
plot(t_cmds, ref(:,1),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('x [m]')
legend('Trilateration','VICON ground truth', 'Command')

subplot(3,1,2)
plot(t_uwb,pos(2,:),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,2),'b','Linewidth',2)
plot(t_cmds, ref(:,2),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('y [m]')
legend('Trilateration','VICON ground truth', 'Command')

subplot(3,1,3)
plot(t_uwb,pos(3,:),'r','Linewidth',2)
grid on
hold on
plot(t_vicon,pos_vicon(:,3),'b','Linewidth',2)
plot(t_cmds, ref(:,3),'--k', 'LineWidth', 1.5)
xlabel('t [s]')
ylabel('z [m]')
legend('Trilateration','VICON ground truth', 'Command')