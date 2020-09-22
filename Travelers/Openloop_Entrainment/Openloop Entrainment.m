clear all
close all
clc
load('Periodic_Solution_JFK_I_1000lux.mat');
Phase_Delay=1;
Periodic_Solution=Periodic_Solution(Phase_Delay*100+1:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
time_shift=[1:23];tol=0.01;
Imax=1000;T=24;D=1600/24;
for i=1:size(time_shift,2)
    n0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),time_shift(i));
    x0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),time_shift(i));
    xc0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),time_shift(i));
    H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),time_shift(i));
    sim('Openloop_Entrainment_I.slx');
    T_openloop(i)=x(end,1);
end
plot(time_shift,T_openloop)
grid on
