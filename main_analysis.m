function z=main_analysis(sd)
 
% ----------------------------------------------------------------------- %
%
% This is the main analyis function, which calls the ODE, fluxes and the
% rest of the functions and plots the results depicted in the main article.
%
% ----------------------------------------------------------------------- %
close all
clear all
%% ------------------ Run Diff_Curien ODEs ------------------------------- %
Y0 = zeros(8,1); %Initial condition
step = 1e-2;
t=0:step:1.5e3;
[T,X] = ode15s(@(t,y) Diff_Curien_sz(t,y,0,0),t,Y0);
 
figure(1)
plot(T,X,'linewidth',2.5)
legend('X_1','X_2','X_3','X_4','X_5','X_6','X_7','X_8')
xlabel('time [min]','FontSize',16)
 
%% ----------- Run Fluxes_Curien to compare with analysis output---------- %
[fluxes,Xdot] = Fluxes_Curien(T,X');
 
figure(2)
plot(T,fluxes,'linewidth',2.5)
legend('v_1','v_2','v_3','v_4','v_5','v_6','v_7','v_8','v_9','v_{10}')
xlabel('time [min]','FontSize',16)
 
%% ------------ Translate the stoichiometry and pathway structure into the
% Stoi matrix below ----------------------------------------------------- %
Stoi = [1 , -1 ,  0 ,  0 ,  0 ,  0 ,  0  ,  0 ,  0 ,  0 
        0 , 1  , -1 , -1 ,  0 ,  0 ,  0  ,  0 ,  0 ,  0
        0 , 0  ,  0 ,  1 , -1 ,  0 ,  0  ,  0 ,  0 ,  0
        0 , 0  ,  1 ,  0 ,  0 , -1 ,  0  ,  0 ,  0 ,  0
        0 , 0  ,  0 ,  0 ,  0 ,  1 , -1  ,  0 ,  0 ,  0
        0 , 0  ,  0 ,  0 ,  0 ,  0 ,  1  , -1 , -1 ,  0
        0 , 0  ,  0 ,  0 ,  0 ,  0 ,  0  ,  0 ,  1 , -1
        0 , 0  ,  0 ,  0 ,  0 ,  0 ,  0  ,  1 ,  0 ,  0];
A = Stoi;    
save Curien_Xdot Xdot A
% ----------- NullSpace is a matrix with each columns constituting a basis 
% for the null sapce, here these are two 10-dimensional vectors, i.e. a 10 
% by 2 NullSpace ------------------------------------------------------- % 
NullSpace = null(Stoi);
 
% NullSpace =
% 
%     0.5374    0.0534
%     0.5374    0.0534
%     0.1162    0.3914
%     0.4212   -0.3380
%     0.4212   -0.3380
%     0.1162    0.3914
%     0.1162    0.3914
%    -0.0000    0.0000
%     0.1162    0.3914
%     0.1162    0.3914
 
% ----------- Gammas are simply the projection of the fluxes on the
% nullspace ------------------------------------------------------------ %
Gammas = NullSpace'*fluxes; 
 
save Curien_data T fluxes Xdot X Gammas A
 
figure(3)
plot(Gammas(1,:),Gammas(2,:),'.')
arrow('start',Gammas(:,1:900:5000)','stop',Gammas(:,2:900:5000)')
 xlabel('\gamma_1','FontSize',16)
 ylabel('\gamma_2','FontSize',16)
 
%%
b = Xdot; % Each row is Xi_dot
figure(4)
plot(Gammas(1,:),Gammas(2,:),'.')
arrow('start',Gammas(:,1:900:5000)','stop',Gammas(:,2:900:5000)')
for j = 1:7
hold on
lb = [-1,-1];
ub = [10.5,Inf];
c = -pinv(A)*b(:,j);
% c = -A\b(:,j);
plotregion(NullSpace,c,lb,ub,'r',0.1)
end
axis equal
xlim([0.5,4.1])
ylim([-0.7,2])
xlabel('\gamma_1','FontSize',16)
ylabel('\gamma_2','FontSize',16)
  
 
%% ------- Finding corners of the non-negative fluxes Gamma-space -------- %
% For numerical reasons only and only for this section of finding the 
% corners, we remove the effect of the fully determined v10 and redefine the
% vectors b consisting of time derivatives of all metabolites but X8 and 
% the 7 by 9 Stoichiometric matrix A9.
 
% The numerical issues are caused because NullSpace has zero entries for v8
% but MATLAB calculation error causes the entries to be a very small on the
% order of 1e-15 but non-zero entries, which introduces an unwanted new
% inequality constraint. 
A9 = [1 , -1 , 0 , 0 , 0 , 0 , 0  , 0 , 0 
      0 , 1  ,-1 ,-1 , 0 , 0 , 0  , 0 , 0 
      0 , 0  , 0 , 1 ,-1 , 0 , 0  , 0 , 0 
      0 , 0  , 1 , 0 , 0 ,-1 , 0  , 0 , 0 
      0 , 0  , 0 , 0 , 0 , 1 ,-1  , 0 , 0 
      0 , 0  , 0 , 0 , 0 , 0 , 1  ,-1 , 0 
      0 , 0  , 0 , 0 , 0 , 0 , 0  , 1 , -1 ];
 
NullSpace9 = null(A9);%(9x2)
v8 = fluxes(8,:);
a = zeros(size(Xdot)); a(6,:) = v8;
b9 = Xdot(1:7,:) + a(1:7,:); 
fluxes9 = fluxes([1:7,9:10],:);
Corners_init = zeros(size(b9,2),2);
j=0;
%%
size(b9)
for k=1:size(b9,2)
    j=j+1;
    cc = [pinv(A9)*b9(:,k);10000];
    AA = [-NullSpace9;1 0];
    [VV,nr,nre]=lcon2vert(AA,cc,[],[]);
    Corners_init(j,:) = VV(1,:);
end
 
figure(4)
plot(Corners_init(:,1),Corners_init(:,2),'.k')
xlabel('\gamma_1','FontSize',16)
ylabel('\gamma_2','FontSize',16)
 
save Corners_init Corners_init
 
%% Adding a constant amount to the fluxes in Set 1 and 2 of fluxes %% 
% Here, two parametrs s and z are introduced to show
% that addition of an equal amount (=s) to certain fluxes belonging to Set 1
% and/or an equal amount (=z) to Set 2 does not change the resulting
% metabolite concentrations. Figure 5 shows how Figure 7A of the main 
% article is plotted. 
s_list=[0,2,4,6,8];
z_list=[0,2,4,6,8];
l=0;
Gam = zeros(2,5,5);
figure(5)
for i=s_list(1:5)
    l=l+1;
    k=0;
    for j=z_list(1:5)
        k=k+1;
        [T,X] = ode15s(@(t,y) Diff_Curien_sz(t,y,i,j),t,Y0);
        fluxes = Fluxes_Curien_sz(T,X',i,j);
        Gammas = NullSpace'*fluxes;
        Gam(:,l,k)= Gammas(:,2500);
        hold on
        ax = gca;
        ax.ColorOrderIndex = 1;
        plot(Gammas(1,:),Gammas(2,:),'.')
    end
end
Gam1 = reshape(Gam(:,1,:),[2,5]);
Gam2 = reshape(Gam(:,:,1),[2,5]);
hold on
plot(Gam1(1,:),Gam1(2,:),'-oc','linewidth',2)
plot(Gam2(1,:),Gam2(2,:),'-or','linewidth',2)
 xlabel('\gamma_1','FontSize',16)
 ylabel('\gamma_2','FontSize',16)
 hold off 
 
%% 
 
end
