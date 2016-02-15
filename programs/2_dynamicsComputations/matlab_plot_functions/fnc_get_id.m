function [id_res] = fnc_get_id (file_id, l_step, r_step, model, id_res, patientWeight)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

id = dlmread(file_id);
id = id/patientWeight;
l_stepDuration = l_step(2)-l_step(1);
r_stepDuration = r_step(2)-r_step(1);

if (~isnan(r_stepDuration))
    hip_r = id(r_step(1):r_step(2),8:10);
    knee_r = id(r_step(1):r_step(2),11:12);
    ankle_r = id(r_step(1):r_step(2),13:14);
    
    hip_r_res = resample(hip_r,100,length(hip_r));
    knee_r_res = resample(knee_r,100,length(knee_r));
    ankle_r_res = resample(ankle_r,100,length(ankle_r));
    
    id_res(model).hip_r_res = hip_r_res;
    id_res(model).knee_r_res = knee_r_res;
    id_res(model).ankle_r_res = ankle_r_res;
else
    id_res(model).hip_r_res = NaN(100,3);
    id_res(model).knee_r_res = NaN(100,2);
    id_res(model).ankle_r_res = NaN(100,2);
end
if (~isnan(l_stepDuration))
    hip_l = id(l_step(1):l_step(2),15:17);
    knee_l = id(l_step(1):l_step(2),18:19);
    ankle_l = id(l_step(1):l_step(2),20:21);
    
    hip_l_res = resample(hip_l,100,length(hip_l));
    knee_l_res = resample(knee_l,100,length(knee_l));
    ankle_l_res = resample(ankle_l,100,length(ankle_l));
    
    id_res(model).hip_l_res = hip_l_res;
    id_res(model).knee_l_res = knee_l_res;
    id_res(model).ankle_l_res = ankle_l_res;
else
    id_res(model).hip_l_res = NaN(100,3);
    id_res(model).knee_l_res = NaN(100,2);
    id_res(model).ankle_l_res = NaN(100,2);
end