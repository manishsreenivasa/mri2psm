function [M,COM,Inertia_com] = fnc_comp_Mass_Com_Inertia(limb, rho, Intervals, Interpolationtype, nSlices) 
% MRI2PSM v0.1 - Patient-specific rigid body modeling from MRI images
% Copyright (c) 2016 Manish Sreenivasa <manish.sreenivasa@iwr.uni-heidelberg.de>, Daniel Gonzalez Alvarado, Carlos Javier Gonzalez Chamorro, Heidelberg University, Germany
%
% Licensed under the zlib license. See LICENSE for more details.
%
% Program: fnc_compute_segment_properties
% The order of the rows is x,y,z,x2,xy,y2,xz,yz,z2

%Obtaining Medium values per slice
B_med = zeros(9,nSlices);
M_med = zeros(9,nSlices);
F_med = zeros(9,nSlices);
for k=1:nSlices;
    if(isempty(limb(k).bonePos)*-1+1)
       N=length(limb(k).bonePos(:,1));
       Xk = limb(k).bonePos(:,1);
       Yk = limb(k).bonePos(:,2);
       Zk = limb(k).bonePos(:,3);
       B_med(1,k) = sum(Xk)/N;
       B_med(2,k) = sum(Yk)/N;         
       B_med(3,k) = sum(Zk)/N;  
       B_med(4,k) = transpose(Xk)*(Xk)/N;
       B_med(5,k) = transpose(Xk)*(Yk)/N;
       B_med(6,k) = transpose(Yk)*(Yk)/N;        
       B_med(7,k) = transpose(Xk)*(Zk)/N;
       B_med(8,k) = transpose(Yk)*(Zk)/N;
       B_med(9,k) = transpose(Zk)*(Zk)/N;   
    end
end
for k=1:nSlices;
    if(isempty(limb(k).musclePos)*-1+1)
       N=length(limb(k).musclePos(:,1));
       Xk = limb(k).musclePos(:,1);
       Yk = limb(k).musclePos(:,2);       
       Zk = limb(k).musclePos(:,3);
       M_med(1,k) = sum(Xk)/N;
       M_med(2,k) = sum(Yk)/N;               
       M_med(3,k) = sum(Zk)/N;  
       M_med(4,k) = transpose(Xk)*(Xk)/N;
       M_med(5,k) = transpose(Xk)*(Yk)/N;
       M_med(6,k) = transpose(Yk)*(Yk)/N;     
       M_med(7,k) = transpose(Xk)*(Zk)/N;
       M_med(8,k) = transpose(Yk)*(Zk)/N;
       M_med(9,k) = transpose(Zk)*(Zk)/N; 
       end
end
for k=1:nSlices;
    if(isempty(limb(k).fatPos)*-1+1)
       N=length(limb(k).fatPos(:,1));
       Xk = limb(k).fatPos(:,1);
       Yk = limb(k).fatPos(:,2);       
       Zk = limb(k).fatPos(:,3);
       F_med(1,k) = sum(Xk)/N;
       F_med(2,k) = sum(Yk)/N;          
       F_med(3,k) = sum(Zk)/N;         
       F_med(4,k) = transpose(Xk)*(Xk)/N;
       F_med(5,k) = transpose(Xk)*(Yk)/N;
       F_med(6,k) = transpose(Yk)*(Yk)/N; 
       F_med(7,k) = transpose(Xk)*(Zk)/N;
       F_med(8,k) = transpose(Yk)*(Zk)/N;
       F_med(9,k) = transpose(Zk)*(Zk)/N; 
    end
end




%Interpolating
Z_itp=linspace(min([limb().sliceZ]),max([limb().sliceZ]),Intervals);
B_med_itp = zeros(9,Intervals);
M_med_itp = zeros(9,Intervals);
F_med_itp = zeros(9,Intervals);
for n=1:9
    B_med_itp(n,:) = interp1([limb().sliceZ], B_med(n,:), Z_itp,Interpolationtype);
    M_med_itp(n,:) = interp1([limb().sliceZ], M_med(n,:), Z_itp,Interpolationtype);
    F_med_itp(n,:) = interp1([limb().sliceZ], F_med(n,:), Z_itp,Interpolationtype);
end
B_area_itp = interp1([limb().sliceZ], [limb().boneArea]  , Z_itp,Interpolationtype);
M_area_itp = interp1([limb().sliceZ], [limb().muscleArea], Z_itp,Interpolationtype);
F_area_itp = interp1([limb().sliceZ], [limb().fatArea]   , Z_itp,Interpolationtype);

%Integrating
dz = (max([limb().sliceZ])-min([limb().sliceZ]))/Intervals;
B_Vintegrals = dz*B_med_itp*transpose(B_area_itp);
M_Vintegrals = dz*M_med_itp*transpose(M_area_itp);
F_Vintegrals = dz*F_med_itp*transpose(F_area_itp);
IntegralsOverMass = [B_Vintegrals, M_Vintegrals, F_Vintegrals]*transpose(rho);

%Obtaining Mass, COM and Inertia
M = dz*[sum( B_area_itp), sum( M_area_itp), sum( F_area_itp)]*transpose(rho);
COM = IntegralsOverMass(1:3)/M;
XXM = IntegralsOverMass(4);
XYM = IntegralsOverMass(5);
YYM = IntegralsOverMass(6);
XZM = IntegralsOverMass(7);
YZM = IntegralsOverMass(8);
ZZM = IntegralsOverMass(9);
Inertia = zeros(3,3);
Inertia(1,:) = [ZZM+YYM, -XYM   , -XZM   ];
Inertia(2,:) = [-XYM   , XXM+ZZM, -YZM   ];
Inertia(3,:) = [-XZM   , -YZM   , XXM+YYM];

%Translating Inertia to COM
XXMcom = XXM-M*(COM(1))*(COM(1));
XYMcom = XYM-M*(COM(1))*(COM(2));
YYMcom = YYM-M*(COM(2))*(COM(2));
XZMcom = XZM-M*(COM(1))*(COM(3));
YZMcom = YZM-M*(COM(2))*(COM(3));
ZZMcom = ZZM-M*(COM(3))*(COM(3));
Inertia_com = zeros(3,3);
Inertia_com(1,:) = [ZZMcom+YYMcom, -XYMcom      , -XZMcom      ];
Inertia_com(2,:) = [-XYMcom      , XXMcom+ZZMcom, -YZMcom      ];
Inertia_com(3,:) = [-XZMcom      , -YZMcom      , XXMcom+YYMcom];

end
