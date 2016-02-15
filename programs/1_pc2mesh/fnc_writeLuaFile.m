function res = fnc_writeLuaFile (luaFileName, model)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

[lua_file] = fopen(luaFileName,'wt');
fprintf(lua_file,'return {\n\tgravity = { 0, 0, -9.81,},\n\tconfiguration = {\n\t\taxis_right = { 0, -1, 0,},\n\t\taxis_front = { 1, 0, 0,},\n\t\taxis_up = { 0, 0, 1,},},\nframes = {\n');
for i=1:length(model)
    fnc_write_lua_body(lua_file, model(i).name, model(i).parent, model(i).rel_transformation, model(i).rel_joint ,...
        model(i).inertia, model(i).com, model(i).mass, model(i).segmentLength, model(i).type, model(i).joint_type, ... ,...
        model(i).marker_names, model(i).marker_values, model(i).mesh_obj, model(i).mesh_dimension, model(i).mesh_center, model(i).mesh_color);
    
end
fprintf(lua_file,'},}');
res = fclose(lua_file);
