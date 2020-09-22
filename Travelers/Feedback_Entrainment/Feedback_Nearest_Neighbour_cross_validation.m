clear all;
clc;
bdclose('all');
global Sleep u H x xc T_opt n I
load('Periodic_Solution_JFK_I_1000lux.mat');
load('Training_Data_10000lux.mat');
Imax=10000;umax=0.2208;
tol=0.01;
Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
Initial_Time=9;% initial local time of subjects
x0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),Initial_Time-6);
xc0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),Initial_Time-6);
H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),Initial_Time-6);
n0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),Initial_Time-6);
Sleep0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,6),Initial_Time-6);
time_shift=[23:-1:1];
for nn=1:size(time_shift,2)
    load('Periodic_Solution_JFK_I_1000lux.mat');
    Time_shift=time_shift(nn);% time shift for international travelers
    Initial_Reference=mod(Initial_Time+Time_shift-6,24);
    [M,N]=min(abs([Periodic_Solution(:,1)-Initial_Reference]));
    Periodic_Solution=Periodic_Solution(N:end,:);
    Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
    sim('Feedback_NearestNeighbor.slx');
    T_Feedback(nn)=x(end,1);
    figure (nn)
    plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'b','linewidth',2)
    hold on
    plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',2)
    hold on
    plot(x(:,1),x(:,2),'k','linewidth',2)
    hold on
    plot(xc(:,1),xc(:,2),'k--','linewidth',2)
    hold on
    plot(u(:,1),u(:,2),'r','linewidth',2)
    hold on
    plot(Sleep(:,1),Sleep(:,2),'g','linewidth',2)
    grid on
    axis([0 u(end,1) -1.3 1.3])
    filename=strcat(num2str(nn),'_Controllable_9am_10000lux_Feedback','.mat');
    save(filename,'Sleep','x','xc','I','u','H','n','Schedule')
end
figure (24)
hold on
time_shift=[1:23];
T_Feedback=[time_shift;T_Feedback];
hold on
plot(T_Feedback(1,:),T_Feedback(2,:))
save('Entrainment_Time_Feedback_9am_10000lux.mat','T_Feedback')

for i=1:23
    if T_opt(i)>T_Feedback(2,i)
        T_opt(i)=T_Feedback(2,i);
    end
end

