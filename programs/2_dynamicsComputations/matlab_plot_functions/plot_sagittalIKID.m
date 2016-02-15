% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Read IK and ID results and plot values

clear; 
figure ('name','Sagittal IK and ID results', 'position', [0 0 1200 500]); 

patientRootDir = '../../../';
patient_info = tdfread([patientRootDir,'sampleData/patient_characteristics.txt'], ',');

r_color = 'b'; l_color = 'r';
width = 1; mnWidth = 2;
[b,a] = butter(3, 0.3, 'low');
fntSz = 10;

modelNames = {'gsm', 'psm'};
lineTypeName = {'--', '-'};
ankleFlexion_offset = [0.5858 0.5338]*180/pi; % Get these values from output of file programs/1_pc2mesh/prog_2_compute_joint_axes

% Load stepping events from corresponding c3d file      
c3d_fileName = [patientRootDir, 'sampleData/c3d/sampleData.c3d'];
[l_step,r_step] = fnc_getSteppingEvents(c3d_fileName);
ik_res = []; id_res = [];

for model = 1:length(modelNames)
    modelType = char(modelNames(model));
    lnSty = char(lineTypeName(model));
    
    ik_fileName = [patientRootDir, 'sampleResults/ik_', modelType, '_sampleData.csv'];
    ik_res = fnc_get_ik (ik_fileName, l_step, r_step, model, ik_res);
        
    id_fileName = [patientRootDir, 'sampleResults/id_', modelType, '_sampleData.txt'];
    id_res = fnc_get_id (id_fileName, l_step, r_step, model, id_res, patient_info.weight(1));
end

timeStamps = [1:100];

ik_res(2).ankle_l_res(:,1) = ik_res(2).ankle_l_res(:,1) + ankleFlexion_offset(1);
ik_res(2).ankle_r_res(:,1) = ik_res(2).ankle_r_res(:,1) + ankleFlexion_offset(2);

subplot(2,6,1); hold on;
plot(timeStamps, -ik_res(1).hip_l_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle','--');
plot(timeStamps, -ik_res(2).hip_l_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle','-');
ylim([-35 50]); title('Left'); ylabel('<-- Ext.       (\circ)       Flex. -->');
text(70, 75, 'Hip Flexion-Extension', 'fontsize', fntSz);
legend('GSM','PSM');

subplot(2,6,2); hold on;
plot(timeStamps, -ik_res(1).hip_r_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, -ik_res(2).hip_r_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle', '-');
ylim([-35 50]); title('Right');
 
subplot(2,6,3); hold on;
plot(timeStamps, ik_res(1).knee_l_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, ik_res(2).knee_l_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle', '-');
ylim([-10 65]); title('Left');
text(70, 88, 'Knee Flexion-Extension', 'fontsize', fntSz);

subplot(2,6,4); hold on;
plot(timeStamps, ik_res(1).knee_r_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, ik_res(2).knee_r_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle', '-');
ylim([-15 65]); title('Right');

subplot(2,6,5); hold on;
plot(timeStamps, -ik_res(1).ankle_l_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, -ik_res(2).ankle_l_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle', '-');
ylim([-25 25]); title('Left');
text(70, 40, 'Ankle Plantar-Dorsiflexion', 'fontsize', fntSz);

subplot(2,6,6); hold on;
plot(timeStamps, -ik_res(1).ankle_r_res(:,1), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, -ik_res(2).ankle_r_res(:,1), 'r', 'linewidth', mnWidth, 'linestyle', '-');
ylim([-25 25]); title('Right');

subplot(2,6,7); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).hip_l_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle','--');
plot(timeStamps, filtfilt(b,a,id_res(2).hip_l_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle','-');
set(gca,'ytick',[-1 -0.5 0 0.5 1]);
ylim([-0.85 0.65]); ylabel('<-- Ext.  (Nm/Kg)   Flex. -->'); xlabel('% Step');

subplot(2,6,8); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).hip_r_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, filtfilt(b,a,id_res(2).hip_r_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle', '-');
set(gca, 'ytick', [-1 -0.5 0 0.5 1]);
ylim([-0.85 0.65]);

subplot(2,6,9); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).knee_l_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, filtfilt(b,a,id_res(2).knee_l_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle', '-');
set(gca, 'ytick', [-1 -0.5 0 0.5 1]);
ylim([-0.5 0.5]); 
 
subplot(2,6,10); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).knee_r_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, filtfilt(b,a,id_res(2).knee_r_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle', '-');
set(gca, 'ytick', [-1 -0.5 0 0.5 1]);
ylim([-0.5 0.5]);

subplot(2,6,11); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).ankle_l_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, filtfilt(b,a,id_res(2).ankle_l_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle', '-');
set(gca, 'ytick', [0 0.5 1 1.5]);
ylim([-0.25 1.25]);
 
subplot(2,6,12); hold on;
plot(timeStamps, filtfilt(b,a,id_res(1).ankle_r_res(:,1)), 'k', 'linewidth', mnWidth, 'linestyle', '--');
plot(timeStamps, filtfilt(b,a,id_res(2).ankle_r_res(:,1)), 'r', 'linewidth', mnWidth, 'linestyle', '-');
set(gca, 'ytick', [0 0.5 1 1.5]);
ylim([-0.25 1.25]);

for i=1:12
    subplot(2,6,i); hold on; set(gca, 'fontsize', fntSz, 'xtick', [0 60 100]);
    plot([0 100],[0 0],'--k');
    plot([60 60],[-100 100],'--k');
    xlim([0 100]);
    axis square;
 end
