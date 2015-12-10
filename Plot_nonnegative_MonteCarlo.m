load NonNegative_V % Contains Y1 Y2 BB VV : 2082 working trajectories in 
% Gamma plane that result in non-negative fluxes
% BB (4x2082) , VV(18738x5001) , Y1,Y2 (5001x2082)
 
load Curien_data_150 % Contains Curien fluxes v1,...,v7,v9,v10 and the 
% corresponding T = 0:0.1:150
% It also has Xdot (7x1501) for t=0:0.1:150 stored in it, as 
% well as matrix A (the 7x9 version)
ind=[1:7,9:10]; %v8 is not included since it's fully determined.
 
% ------------------------------ A number of flux sets that are reletively 
% visually diverse are selected for plotting ---------------------------- %
smallsetind = [71,619,1417,49,12,45,11,2];
 
figure(6)
for i=1:9
subplot(3,3,i)
plot(T,VV(i:9:end,:))
hold on
plot(T,fluxes(i,:)+VV(i,1)-fluxes(i,1),'.k')
xlim([0,150])
str = sprintf('v_{%d}',ind(i));
if (i>6)&&(i<10)
xlabel('time [$s$]','Interpreter','latex')
end
if (i==1)|(i==4)|(i==7)
ylabel('Flux [$\mu M^{-1}.s^{-1}$]','Interpreter','latex')
end
title(str);
hold off
end
