function fnc_plot3LineFromTo(from, to, col, wd)
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%

if nargin < 4
    wd=1;
end
if nargin < 3
    col='k';
end
[a,b]=size(from);
for i=1:a
    plot3(from(i,1), from(i,2), from(i,3), 'linestyle','none','marker','o','color',col);
    plot3(to(i,1), to(i,2), to(i,3), 'linestyle','none','marker','o','color',col);
    plot3([from(i,1) to(i,1)],[from(i,2) to(i,2)], [from(i,3) to(i,3)], 'linestyle', '-', 'color', col,'linewidth',wd);
end
