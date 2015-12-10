load Curien_data % Contains T fluxes Xdot X Gammas A
load Corners_init % Contains Corners_init
%% ------------ Quadratic programing: minimizing L2 norm --------------- %%
v_quad = zeros(10,size(Xdot,2));
H=eye(10);
f=zeros(10,1);
for i=1:size(Xdot,2) 
v_quad(:,i) = quadprog(H,f,[],[],A,Xdot(:,i),zeros(10,1),[]);
end
Gammas_quad = NullSpace'*v_quad; 
 
save min_norm_solution v_quad Gammas_quad T
%% ---- Quadratic prog. only for a subset of time points -> faster ---- %%
v_quad2 = zeros(10,size(Xdot,2));
H=eye(10);
f=zeros(10,1);
for i=1:100:size(Xdot,2) 
v_quad2(:,i) = quadprog(H,f,[],[],A,Xdot(:,i),zeros(10,1),[]);
end
V0_2 = v_quad2(:,1:100:size(Xdot,2))-pinv(A)*Xdot(:,1:100:size(Xdot,2)); 
Gammas_quad2 = NullSpace'*V0_2; 
%% ------------ Linear programing: minimizing L1 norm --------------- %%
v_lin = zeros(10,size(Xdot,2));
f=ones(10,1);
for i=1:100:size(Xdot,2) 
v_lin(:,i) = linprog(f,[],[],A,Xdot(:,i),zeros(10,1),[]);
end
% V0_lin = v_lin(:,1:100:size(Xdot,2))-pinv(A)*Xdot(:,1:100:size(Xdot,2)); 
Gammas_lin = NullSpace'*v_lin(:,1:100:size(Xdot,2)); 
%%
figure(7)
hold on
plot(Gammas(1,:),Gammas(2,:),'.')
plot(Gammas_quad(1,:),Gammas_quad(2,:),'.')
plot(Gammas_lin(1,:),Gammas_lin(2,:),'.')
plot(Corners_init(:,1),Corners_init(:,2),'.k')
 
legend('Actual Gamma trajectory','Lin. prog','Quad. prog','Cone corners')
xlabel('\gamma_1','FontSize',16)
ylabel('\gamma_2','FontSize',16)
hold off
