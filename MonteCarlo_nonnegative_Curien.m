% For plotting purposes only, the trajectories are shifted such that they all start from Gamma(0)
Y0=[16,3];
    
% This will be the matrix containing the i-th B resulting in 
% non-negative fluxes in the form of the vector [b11,b12,b21,b22]^T as its 
% i-th column. 
Y1 = []; % This will contain the corresponding Gamma1 as its i-th column
Y2 = []; % This will contain the corresponding Gamma2 as its i-th column
VV = []; % And finally the corresponding fluxes vs. time are stored in VV. 
% If this script finds m set of non-negative fluxes, then BB is 4 by m, Y1
% and Y2 are l by m, and VV is an l by m*10 
 
load Curien_data % This has Xdot (7x1501) for T=0:0.1:150 stored in it, as 
% well as matrix A (the 7x9 version)
% It also contains Curien fluxes v1,...,v7,v9,v10, which are not used in
% this script.
 
step = 1e-1;
t=0:step:1.5e2;
NullSpace = null(A);
 
for i = 1:1e5
    b11 = -1+2*rand; % in [-1,1]
    b22 = -b11-rand; % in [1-b11,b11]
    b12 = -1+2*rand; % in [-1,1]
    b21 = -1+2*rand; % in [-1,1]
    % This choice of b22 results in trace(B) = b11+b22 <=0
    deter=b11*b22-b12*b21;
    
    if deter>0 
        B=[b11 b12;b21 b22];
        [T,Y] = ode45(@(t,y) Diff_gamma(t,y,B),t,Y0);
        V = pinv(A)*Xdot + NullSpace*Y';
        if min(min(V))>= -1e-4
            disp('Found one')
            Y1 =[Y1,Y(:,1)];
            Y2 =[Y2,Y(:,2)];
            BB=[BB,reshape(B,4,1)];
            VV=[VV;V];
        end
    end
end
save NonNegative_V Y1 Y2 BB VV T
