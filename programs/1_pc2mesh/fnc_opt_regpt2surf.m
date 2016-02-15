function [C, dist0, flag_constraint_saturation] = fnc_opt_regpt2surf (node, elem, p, maxiter, lambda, trivial_cost, pmask, C0, C_ub, C_lb, cmask, verbose)
%  [A, b, newpos] = fnc_opt_regpt2surf (node, elem, p, maxiter, lambda, trivial_cost, pmask, C0, C_ub, C_lb, cmask, verbose)
%  Perform point cloud registration to a triangular surface
%  (surface can be either triangular or cubic), Gauss-Newton method
%  is used for the calculation
%
%  Original author: Qianqian Fang <fangq at nmr.mgh.harvard.edu>
%  date: 12/12/2008
%  Please find more information at http://iso2mesh.sf.net/cgi-bin/index.cgi?metch
%  The orginal function was part of "metch" toobox, issued under the GNU GPLv2 license http://www.gnu.org/licenses/gpl-2.0.html
%
% Modified 16.01.2016 by Manish Sreenivasa <manish.sreenivasa at iwr.uni-heidelberg.de>
%    - modified function interface
%    - included contraints on affine transformation matrix
%    - introduced homotopy variable while computing residuals
% Modified code released as part of the mri2psm toolbox
% https://github.com/manishsreenivasa/mri2psm
%                                                                    
% parameters:
%      node: node coordinate of the surface mesh (nn x 3)
%      elem: element list of the surface mesh (3 columns for
%            triangular mesh, 4 columns for cubic surface mesh)
%      p: points to be registered, 3 columns for x,y and z respectively
%      maxiter: an integer, specifying the optimization iterations
%      lambda: homotopy variable
%      trivial_cost: default cost for lambda = 0
%      pmask: a mask vector with the same length as p, determines the
%         method to handle the point, if pmask(i)=-1, the point is a free
%         node and can be move by the optimization, if pmask(i)=0, the
%         point is fixed; if pmask(i)=n>0, the distance between p(i,:)
%         and node(n,:) will be part of the object function and be optimized
%      C0: a 12x1 vector [A0(:);b0], as the initial guess for the affine transformation
%      C_ub: a 12x1 vector defining upper bounds of C
%      C_lb: a 12x1 vector defining lower bounds of C
%      cmask: a binary 12x1 vector, determines which element of [A(:);b] will be optimized
%          if cmask(i)=0, the corresponding coefficient will not be updated
%      verbose: print iteration results if 1
%
% outputs:
%      C: 12x1 vector, the updated affine transformation
%      dist0: residual errors
%      flag_constraint_saturation:whether or not iterations were stopped due to hitting limits


C = C0;
C_update = C0;
flag_constraint_saturation = 0;

delta = 1e-4;
newpos = (reshape(C(1:9),3,3)*p'+repmat(C(end-2:end),1,size(p,1)))';
nv = nodesurfnorm(node,elem);
clen = length(C);
cuplist = find(cmask == 1);
pfree = find(pmask < 0);
pfix = find(pmask > 0);

% Start Gauss-Newton iterations
for iter=1:maxiter
    
    % calculate the current residual: the sum of distances to the surface
    dist0 = zeros(length(pfree)+length(pfix),1);        
    if(~isempty(pfree))
        [dist0(pfree),cn0] = dist2surf(node,nv,newpos(pfree,:));
    end
    if(~isempty(pfix))
        fixdist = node(pmask(pfix),:) - newpos(pfix,:);
        dist0(pfix) = sqrt(sum(fixdist.*fixdist,2));
    end
    
    % Compute homotopic cost
    dist0 = (1-lambda).*trivial_cost + lambda.*dist0;
    
    if verbose
        fprintf('iter = %d error = %f\n',iter,sum(abs(dist0)));
    end
    
    % build the Jacobian (sensitivity) matrix
    J = zeros(length(dist0),length(C));
    for i = 1:clen
        if (cmask(i) == 0) continue; end
        dC=C;
        if(C(i))
            dC(i)=C(i)*(1+delta);
        else
            dC(i)=C(i)+delta;
        end
        newp=(reshape(dC(1:9),3,3)*p'+repmat(dC(end-2:end),1,size(p,1)))';
        
        dist=zeros(length(pfree)+length(pfix),1);
        if(~isempty(pfree))
            if(length(cn0)==length(pfree))
                dist(pfree)=dist2surf(node,nv,newp(pfree,:),cn0);
            else
                dist(pfree)=dist2surf(node,nv,newp(pfree,:));
            end
        end
        if(~isempty(pfix))
            fixdist=node(pmask(pfix),:)-newp(pfix,:);
            dist(pfix)=sqrt(sum(fixdist.*fixdist,2));
        end
   
        % J = dL/dC
        J(:,i)=(dist-dist0)/(dC(i)-C(i));
    end
    
    % Weight the matrix (normalization)
    wj = sqrt(sum(J.*J));
    J = J./repmat(wj,length(dist0),1);
    
    % Calculate the update: J*dC = dL
    dC = (J\dist0)./wj';
    C_update(cuplist) = C(cuplist) - 0.5*dC(cuplist);
        
    % If C_update is beyond limits, then break
    if ~isempty(find (C_update > C_ub)) || ~isempty(find (C_update < C_lb))
        flag_constraint_saturation = 1;
        break;
    else
        C(cuplist) = C_update(cuplist);
    end
        
    % Get the updated positions with the calculated A and b
    newpos = (reshape(C(1:9),3,3)*p'+repmat(C(end-2:end),1,size(p,1)))';
end
