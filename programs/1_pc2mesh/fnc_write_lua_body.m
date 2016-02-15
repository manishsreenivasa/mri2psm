function fnc_write_lua_body(fid, Name, Parent, Ematrix ,Joint_position,Imatrix,COMmatrix ,Mass ,Length , Body_name ,Join_type , Marker_names, Marker_positions, Object_file, Dimensions, Mesh_center, Color, bWrite_Visuals)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

if nargin < 18
    bWrite_Visuals = 1;
end

fprintf(fid,'\n{name   = "%s",',Name);
fprintf(fid,'\n\tparent = "%s",',Parent);
fprintf(fid,'\n\tjoint_frame = {');
fprintf(fid,'\n\t\tr = '); fnc_write_lua_matrix(fid,Joint_position ,1);
fprintf(fid,'\n\t\tE = \n'); fnc_write_lua_matrix(fid,Ematrix,3);   fprintf(fid,'},');
fprintf(fid,'\n\tbody = {');
fprintf(fid,'\n\t\tname   = "%s",',Body_name);
fprintf(fid,'\n\t\tmass   = %f,\n\t\tlength = %f,',Mass,Length);
fprintf(fid,'\n\t\tcom = '); fnc_write_lua_matrix(fid,COMmatrix,1);
fprintf(fid,'\n\t\tinertia = \n'); fnc_write_lua_matrix(fid,Imatrix,3); fprintf(fid,'},');
fprintf(fid,'\n\tjoint= \n');
[n1,n2] = size(Join_type);
if n1==1
    fprintf(fid,'{');fnc_write_lua_matrix(fid,Join_type,2); fprintf(fid,'},');
else
    fnc_write_lua_matrix(fid,Join_type,2);
end
if bWrite_Visuals
    fprintf(fid,'\n\tmarkers = {');
    [MarkN,MarkM]=size(Marker_names);
    if size(MarkN)>0
        for l=1:MarkN
            fprintf(fid,'\n\t\t%s = ',Marker_names(l,:)); 
            fnc_write_lua_matrix(fid,Marker_positions(l,:),1)
        end
    end
    fprintf(fid,'},');
    
    fprintf(fid,'\n\tvisuals = {{');
    fprintf(fid,'\n\t\tname        = "%s",',Body_name);
    fprintf(fid,'\n\t\tsrc         = "%s",',Object_file);
    fprintf(fid,'\n\t\tdimensions  = '); fnc_write_lua_matrix(fid,Dimensions,1);
    fprintf(fid,'\n\t\tmesh_center = '); fnc_write_lua_matrix(fid,Mesh_center,1);
    fprintf(fid,'\n\t\tcolor       = '); fnc_write_lua_matrix(fid,Color,1); fprintf(fid,'  },},},\n');
else
    fprintf(fid,' },\n');
end
