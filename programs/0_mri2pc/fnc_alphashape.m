function [V,SHP] = fnc_alphashape(X,R,fig)
%   Alpha shape of 2D or 3D point set.
%   V = fnc_alphavol(X,R) gives the area or volume V of the basic alpha shape
%   for a 2D or 3D point set. X is a coordinate matrix of size Nx2 or Nx3.
%
%   R is the probe radius with default value R = Inf. In the default case
%   the basic alpha shape (or alpha hull) is the convex hull.
%
%   [V,S] = fnc_alphavol(X,R) outputs a structure S with fields:
%    S.tri - Triangulation of the alpha shape (Mx3 or Mx4)
%    S.vol - Area or volume of simplices in triangulation (Mx1)
%    S.rcc - Circumradius of simplices in triangulation (Mx1)
%    S.bnd - Boundary facets (Px2 or Px3)
%
%   fnc_alphavol(X,R,1) plots the alpha shape.
%   Author: Jonas Lundgren <splinefit@gmail.com> 2010
%
% Copyright (c) 2010, Jonas Lundgren
% All rights reserved.
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%    * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%      * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% Modified 16.01.2016, Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>
%	- Changed function name to follow scheme of the MRI2PSM toolbox
%

if nargin < 2 || isempty(R), R = inf; end
if nargin < 3, fig = 0; end

% Check coordinates
dim = size(X,2);
if dim < 2 || dim > 3
    error('alphavol:dimension','X must have 2 or 3 columns.')
end

% Check probe radius
if ~isscalar(R) || ~isreal(R) || isnan(R)
    error('alphavol:radius','R must be a real number.')
end

% Unique points
[X,imap] = unique(X,'rows');

% Delaunay triangulation
T = delaunay(X);

% Remove zero volume tetrahedra since
% these can be of arbitrary large circumradius
if dim == 3
    n = size(T,1);
    vol = volumes(T,X);
    epsvol = 1e-12*sum(vol)/n;
    T = T(vol > epsvol,:);
    holes = size(T,1) < n;
end

% Limit circumradius of simplices
[~,rcc] = circumcenters(TriRep(T,X));
T = T(rcc < R,:);
rcc = rcc(rcc < R);

% Volume/Area of alpha shape
vol = volumes(T,X);
V = sum(vol);

% Return?
if nargout < 2 && ~fig
    return
end

% Turn off TriRep warning
warning('off','MATLAB:TriRep:PtsNotInTriWarnId')

% Alpha shape boundary
if ~isempty(T)
    % Facets referenced by only one simplex
    B = freeBoundary(TriRep(T,X));
    if dim == 3 && holes
        % The removal of zero volume tetrahedra causes false boundary
        % faces in the interior of the volume. Take care of these.
        B = trueboundary(B,X);
    end
else
    B = zeros(0,dim);
end

% Plot alpha shape
if fig
    if dim == 2
        % Plot boundary edges and point set
        x = X(:,1);
        y = X(:,2);
        plot(x(B)',y(B)','r','linewidth',2), hold on
        plot(x,y,'k.'), hold off
        str = 'Area';
    elseif ~isempty(B)
        % Plot boundary faces
        trisurf(TriRep(B,X),'FaceColor','red','FaceAlpha',1/3);
        str = 'Volume';
    else
        cla
        str = 'Volume';
    end
    axis equal
    str = sprintf('Radius = %g,   %s = %g',R,str,V);
    title(str,'fontsize',12)
end

% Turn on TriRep warning
warning('on','MATLAB:TriRep:PtsNotInTriWarnId')

% Return structure
if nargout == 2
    SHP = struct('tri',imap(T),'vol',vol,'rcc',rcc,'bnd',imap(B));
end


%--------------------------------------------------------------------------
function vol = volumes(T,X)
%VOLUMES Volumes/areas of tetrahedra/triangles

% Empty case
if isempty(T)
    vol = zeros(0,1);
    return
end

% Local coordinates
A = X(T(:,1),:);
B = X(T(:,2),:) - A;
C = X(T(:,3),:) - A;
    
if size(X,2) == 3
    % 3D Volume
    D = X(T(:,4),:) - A;
    BxC = cross(B,C,2);
    vol = dot(BxC,D,2);
    vol = abs(vol)/6;
else
    % 2D Area
    vol = B(:,1).*C(:,2) - B(:,2).*C(:,1);
    vol = abs(vol)/2;
end


%--------------------------------------------------------------------------
function B = trueboundary(B,X)
%TRUEBOUNDARY True boundary faces
%   Remove false boundary caused by the removal of zero volume
%   tetrahedra. The input B is the output of TriRep/freeBoundary.

% Surface triangulation
facerep = TriRep(B,X);

% Find edges attached to two coplanar faces
E0 = edges(facerep);
E1 = featureEdges(facerep, 1e-6);
E2 = setdiff(E0,E1,'rows');

% Nothing found
if isempty(E2)
    return
end

% Get face pairs attached to these edges
% The edges connects faces into planar patches
facelist = edgeAttachments(facerep,E2);
pairs = cell2mat(facelist);

% Compute planar patches (connected regions of faces)
n = size(B,1);
C = sparse(pairs(:,1),pairs(:,2),1,n,n);
C = C + C' + speye(n);
[~,p,r] = dmperm(C);

% Count planar patches
iface = diff(r);
num = numel(iface);

% List faces and vertices in patches
facelist = cell(num,1);
vertlist = cell(num,1);
for k = 1:num

    % Neglect singel face patches, they are true boundary
    if iface(k) > 1
        
        % List of faces in patch k
        facelist{k} = p(r(k):r(k+1)-1);
        
        % List of unique vertices in patch k
        vk = B(facelist{k},:);
        vk = sort(vk(:))';
        ik = [true,diff(vk)>0];
        vertlist{k} = vk(ik);
        
    end
end

% Sort patches by number of vertices
ivert = cellfun(@numel,vertlist);
[ivert,isort] = sort(ivert);
facelist = facelist(isort);
vertlist = vertlist(isort);

% Group patches by number of vertices
p = [0;find(diff(ivert));num] + 1;
ipatch = diff(p);

% Initiate true boundary list
ibound = true(n,1);

% Loop over groups
for k = 1:numel(ipatch)

    % Treat groups with at least two patches and four vertices
    if ipatch(k) > 1 && ivert(p(k)) > 3

        % Find double patches (identical vertex rows)
        V = cell2mat(vertlist(p(k):p(k+1)-1));
        [V,isort] = sortrows(V);
        id = ~any(diff(V),2);
        id = [id;0] | [0;id];
        id(isort) = id;

        % Deactivate faces in boundary list
        for j = find(id')
            ibound(facelist{j-1+p(k)}) = 0;
        end
        
    end
end

% Remove false faces
B = B(ibound,:);
