 
load min_norm_solution % Contains v_quad, Gammas_quad, and T
load Curien_data % Contains Curien fluxes, Xdot, X, and Gammas as well aso 
% the corresponding T = 0:0.01:1500, and the matrix A (the 8x10 version)
 
index = [1:7,9:10]; %v8 is not included since it's fully determined.
v_minnorm = v_quad(index,:);
V = fluxes(index,:);
 
%% ------------------------------------------------------------------%
% -------------------- Plotting fluxes vs. time ---------------------%
% -------------------------------------------------------------------%
figure(8)
for i=1:9
subplot(3,3,i)
plot(T,V(i,:))
hold on
plot(T,v_minnorm(i,:))
 
str = sprintf('v_{%d}',ind(i));
if (i>6)&&(i<10)
    xlabel('time [$s$]','Interpreter','latex')
end
if (i==1)|(i==4)|(i==7)
    ylabel('Flux [$\mu M^{-1}.s^{-1}$]','Interpreter','latex')
end
if (i==1)
    legend('Actual','Min-energy')
end
title(str);
hold off
xlim([0 1500])
end
%% ------------------------------------------------------------------%
% ---------------- Plotting fluxes vs. their substrates -------------%
% -------------------------------------------------------------------%
% ---------------- One-substrate fluxes -----------------------------%
figure(9)
subplot(4,1,1)
title('$v_5$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm(5,:))
hold on
plot(X(:,3),V(5,:))
hold off
subplot(4,1,2)
hold on
title('$v_6$ $vs.$ $X_4$','Interpreter','latex')
plot(X(:,4),v_minnorm(6,:))
hold on
plot(X(:,4),V(6,:))
hold off
subplot(4,1,3)
hold on
title('$v_7$ $vs.$ $X_5$','Interpreter','latex')
plot(X(:,5),v_minnorm(7,:))
hold on
plot(X(:,5),V(7,:))
hold off
subplot(4,1,4)
hold on
title('$v_{10}$ $vs.$ $X_7$','Interpreter','latex')
plot(X(:,7),v_minnorm(9,:))
hold on
plot(X(:,7),V(9,:))
hold off
 
%% ------------------ Two-substrate fluxes --------------------------- %
figure(10)
subplot(4,1,1)
title('$v_5$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm_new(5,:))
hold on
plot(X(:,3),V(5,:))
hold off
subplot(4,1,2)
hold on
title('$v_6$ $vs.$ $X_4$','Interpreter','latex')
plot(X(:,4),v_minnorm_new(6,:))
hold on
plot(X(:,4),V(6,:))
hold off
subplot(4,1,3)
hold on
title('$v_7$ $vs.$ $X_5$','Interpreter','latex')
plot(X(:,5),v_minnorm_new(7,:))
hold on
plot(X(:,5),V(7,:))
hold off
subplot(4,1,4)
hold on
title('$v_{10}$ $vs.$ $X_7$','Interpreter','latex')
plot(X(:,7),v_minnorm_new(9,:))
hold on
plot(X(:,7),V(9,:))
hold off
%% ---------------- Three-substrate flux v1 ------------------------- %
figure(11)
subplot(1,3,1)
title('$v_1$ $vs.$ $X_1$','Interpreter','latex')
hold on
plot(X(:,1),v_minnorm(1,:))
hold on
plot(X(:,1),V(1,:))
hold off
subplot(1,3,2)
title('$v_1$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm(1,:))
hold on
plot(X(:,3),V(1,:))
hold off
subplot(1,3,3)
title('$v_1$ $vs.$ $X_6$','Interpreter','latex')
hold on
plot(X(:,6),v_minnorm(1,:))
hold on
plot(X(:,6),V(1,:))
hold off
f=gcf;
rez=300;
set(f,'paperunits','inches','paperposition',[0 0 4.5 1]); %dont need to change anything here
path= 'C:\Users\sdolatshahi3.AD\Dropbox\0- 3rd journal-identifiability\PaperWriteup\MATLAb files'; %the folder where you want to put the file
name='plot15.png'; %what you want the file to be called
print(f,fullfile(path,name),'-dpng',['-r',num2str(rez)],'-opengl') %save fil
 
%% ------------------------------------------------------------------%
% ---------- Removing the folding-over phenomenon for v6  ---------- % 
% -------------------------------------------------------------------%
[Y,Ind] = max(v_minnorm(6,:));
X4 = X(:,4);
v6 = v_minnorm(6,1:Ind);
v6q = interp1(X4(1:Ind)',v6,X4(Ind+1:end),'linear','extrap'); 
 
figure(12)
subplot(1,2,1)
plot(X4(1:Ind)',v6,'.')
hold on
plot(X4(Ind+1:end)',v_minnorm(6,Ind+1:end),'.r')
hold on
plot(X4(Ind+1:end),v6q,'.g')
legend('v_6','folded-over section of v_6','mapped to initial branch','location','northwest')
xlabel('$X_4$ [$\mu M^{-1}$]','Interpreter','latex')
ylabel('$v_6$ flux [$\mu M^{-1}.s^{-1}$]','Interpreter','latex')
hold  off
 
subplot(1,2,2)
plot(T,V(6,:),'.r')
hold on
plot(T,[v6,v6q'],'.')
legend('v_6',' v_{6.new}')
xlabel('time [$s$]','Interpreter','latex')
ylabel('$v_6$ flux [$\mu M^{-1}.s^{-1}$]','Interpreter','latex')
hold off
 
v6_new = [v6,v6q'];
 
Xdot_new = Xdot+[0;0;0;1;-1;0;0;0]*v6_new; 
A_new = [A(:,1:5),A(:,7:10)];
V_new = [V(1:5,:);V(7:9,:)];
size(Xdot_new)
size(A_new)
size(V_new)
 
%% -------------------------------------------------------------------%%
% -----------------  Re-calculating v_minnorm ------------------------ % 
% -------------------------------------------------------------------  %
v_quad_new = zeros(9,size(Xdot_new,2));
H=eye(9);
f=zeros(9,1);
for i=1:size(Xdot_new,2) 
v_quad_new(:,i) = quadprog(H,f,[],[],A_new,Xdot_new(:,i),-1e-7*ones(9,1),[]); %-1e-7 is set to overcome  
end
Gammas_quad_new = NullSpace'*v_quad_new; 
 
%% ------------------------------------------------------------------%
% -------------- The corrected fluxes vs. time --------------------- %
% -------------------------------------------------------------------%
v_temp = [v_quad_new(1:5,:);v6_new;v_quad_new(6:9,:)];
v_minnorm_new = v_temp(index,:);
 
figure(13)
for i=1:9
subplot(3,3,i)
plot(T,V(i,:))
hold on
plot(T,v_minnorm_new(i,:))
 
str = sprintf('v_{%d}',ind(i));
if (i>6)&&(i<10)
    xlabel('time [$s$]','Interpreter','latex')
end
if (i==1)|(i==4)|(i==7)
    ylabel('Flux [$\mu M^{-1}.s^{-1}$]','Interpreter','latex')
end
if (i==1)
    legend('Actual','Min-energy')
end
title(str);
hold off
xlim([0 1500])
end
 
%% ------------------------------------------------------------------%
% ---------- Plotting corrected fluxes vs. their substrates ---------%
% -------------------------------------------------------------------%
% ---------------- One-substrate fluxes -----------------------------%
figure(14)
subplot(4,1,1)
title('$v_5$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm_new(5,:))
hold on
plot(X(:,3),V(5,:))
hold off
subplot(4,1,2)
hold on
title('$v_6$ $vs.$ $X_4$','Interpreter','latex')
plot(X(:,4),v_minnorm_new(6,:))
hold on
plot(X(:,4),V(6,:))
hold off
subplot(4,1,3)
hold on
title('$v_7$ $vs.$ $X_5$','Interpreter','latex')
plot(X(:,5),v_minnorm_new(7,:))
hold on
plot(X(:,5),V(7,:))
hold off
subplot(4,1,4)
hold on
title('$v_{10}$ $vs.$ $X_7$','Interpreter','latex')
plot(X(:,7),v_minnorm_new(9,:))
hold on
plot(X(:,7),V(9,:))
hold off
 
%% ------------------ Two-substrate fluxes --------------------------- %
figure(15)
subplot(4,2,1)
title('$v_2$ $vs.$ $X_1$','Interpreter','latex')
hold on
plot(X(:,1),v_minnorm_new(2,:))
hold on
plot(X(:,1),V(2,:))
hold off
subplot(4,2,2)
title('$v_2$ $vs.$ $X_2$','Interpreter','latex')
hold on
plot(X(:,2),v_minnorm_new(2,:))
hold on
plot(X(:,2),V(2,:))
hold off
subplot(4,2,3)
title('$v_3$ $vs.$ $X_2$','Interpreter','latex')
hold on
plot(X(:,2),v_minnorm_new(3,:))
hold on
plot(X(:,2),V(3,:))
hold off
subplot(4,2,4)
title('$v_3$ $vs.$ $X_6$','Interpreter','latex')
hold on
plot(X(:,6),v_minnorm_new(3,:))
hold on
plot(X(:,6),V(3,:))
hold off
subplot(4,2,5)
title('$v_4$ $vs.$ $X_2$','Interpreter','latex')
hold on
plot(X(:,2),v_minnorm_new(4,:))
hold on
plot(X(:,2),V(4,:))
hold off
subplot(4,2,6)
title('$v_4$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm_new(4,:))
hold on
plot(X(:,3),V(4,:))
hold off
subplot(4,2,7)
title('$v_9$ $vs.$ $X_6$','Interpreter','latex')
hold on
plot(X(:,6),v_minnorm_new(8,:))
hold on
plot(X(:,6),V(8,:))
hold off
subplot(4,2,8)
title('$v_9$ $vs.$ $X_7$','Interpreter','latex')
hold on
plot(X(:,7),v_minnorm_new(8,:))
hold on
plot(X(:,7),V(8,:))
hold off
 
%% --------------- Three-substrate fluxes --------------------------- %
figure(16)
subplot(1,3,1)
title('$v_1$ $vs.$ $X_1$','Interpreter','latex')
hold on
plot(X(:,1),v_minnorm_new(1,:))
hold on
plot(X(:,1),V(1,:))
hold off
subplot(1,3,2)
title('$v_1$ $vs.$ $X_3$','Interpreter','latex')
hold on
plot(X(:,3),v_minnorm_new(1,:))
hold on
plot(X(:,3),V(1,:))
hold off
subplot(1,3,3)
title('$v_1$ $vs.$ $X_6$','Interpreter','latex')
hold on
plot(X(:,6),v_minnorm_new(1,:))
hold on
plot(X(:,6),V(1,:))
hold off

