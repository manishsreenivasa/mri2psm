function []= fnc_write_lua_matrix(file,matrix,tabnumber)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

[N,M]=size(matrix);
if N>0
    for k=1:N
        for n=1:tabnumber
             fprintf(file,'\t');
        end
      if k==1 && N>1
          fprintf(file,'{');
      end
      fprintf(file,'{');
      for j=1:M
      	fprintf(file,'  %f  ,',matrix(k,j));
      end
      if k==N
      	if N>1
            fprintf(file,'},},');
        else
            fprintf(file,'},');
      	end
      else
      	fprintf(file,'},\n');
      end
    end
else
     for n=1:tabnumber
     	fprintf(file,'\t');
     end
    fprintf(file,'{},');
end
end
