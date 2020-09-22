% generate and plot the periodic solutions of the three two-process models 
% periodic solution of full KFK model with sleep
clear all
clc
T=24;
D=1600/24;
Imax=1000;
n0=0;
x0=1.06;xc0=0;H0=0.4;
sim('JFK_Sleep_I.slx');
Periodic_Solution=[x(160801:300001,:),xc(160801:300001,2),H(160801:300001,2),I(160801:300001,2),sleep(160801:300001,2),n(160801:300001,2)];
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
subplot(3,1,1)
plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'b','linewidth',3);
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b-.','linewidth',3);
grid on
axis([0 48  -1.5 1.5])
subplot(3,1,2)
plot(Periodic_Solution(:,1),Periodic_Solution(:,4),'b','linewidth',3);
axis([0 48  0 1])
grid on
subplot(3,1,3)
plot(Periodic_Solution(:,1),mod(atan2(-Periodic_Solution(:,3),Periodic_Solution(:,2)),2*pi),'b','linewidth',3)
axis([0 48  0 8])
grid on
%save('Periodic_Solution_JFK_I_1000lux.mat','Periodic_Solution')

% periodic solution of 2nd-order JFK model with sleep
clear all
clc
T=24;
D=1600/24;
umax=0.1731;
x0=1.06;xc0=0;H0=0.4;
sim('JFK_Sleep_u.slx');
Periodic_Solution=[x(160801:300001,:),xc(160801:300001,2),H(160801:300001,2),u(160801:300001,2),sleep(160801:300001,2)];
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
subplot(3,1,1)
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'r','linewidth',3);
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'r-.','linewidth',3);
grid on
legend('x_{3rd}','x_{c3rd}','x_{2nd}','x_{c2nd}');
subplot(3,1,2)
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,4),'r','linewidth',3);
subplot(3,1,3)
hold on
plot(Periodic_Solution(:,1),mod(atan2(-Periodic_Solution(:,3),Periodic_Solution(:,2)),2*pi),'r','linewidth',3)
%save('Periodic_Solution_JFK_u_1000lux.mat','Periodic_Solution')

% plot the periodic solution of phase reduced model with sleep
clear all
clc
load('f_Kronauer.mat');
T=24;
D=1600/24;
umax=0.1731;
theta0=0;H0=0.4;
sim('PRC_Sleep_u.slx');
Periodic_Solution=[theta(160801:300001,:),H(160801:300001,2),u(160801:300001,2),sleep(160801:300001,2)];
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
subplot(3,1,2)
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'k','linewidth',3);
legend('\theta_{3rd}','\theta_{2nd}','\theta_{1st}');
subplot(3,1,3)
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'k','linewidth',3);
legend('H_{3rd}','H_{2nd}','H_{1st}');
%save('Periodic_Solution_PRC_u_1000lux.mat','Periodic_Solution')
