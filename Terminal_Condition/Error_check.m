clear all
close all
clc
tol=0.01;
T=[];error=[];
for i=1:1000
    N=round(2400*rand);
    if N==0 
        N=1;
    end
    load('Periodic_Solution_JFK_I_1000lux.mat');
    Periodic_Solution=Periodic_Solution(N:end,:);
    Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);    
    Theta=2*pi*rand;
    Phi=pi*rand;
    x0=Periodic_Solution(1,2)+sqrt(tol)*sin(Phi)*cos(Theta);
    xc0=Periodic_Solution(1,3)+sqrt(tol)*sin(Phi)*sin(Theta);
    H0=Periodic_Solution(1,4)+sqrt(tol)*cos(Phi);
    n0=rand;
    Sleep0=Periodic_Solution(1,6);
    sim('Open_Loop_shift_notol_initial.slx');
    T=[T;[max(Error(:,2))]];
    hold on
    plot(Error(:,1),Error(:,2))
    grid on
end
histogram(T);

