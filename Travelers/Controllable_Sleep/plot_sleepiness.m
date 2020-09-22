% plot sleepiness
clear all
close all
clc
time_shift=16;
cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Spontaneous_Sleep'
filename=strcat(num2str(time_shift),'_shift_Spontaneous_Sleep_7am_1000lux','.mat');
load(filename);
plot(x(:,1),H(:,2)-0.1333*x(:,2),'b','linewidth',6)
axis([0 x(end,1) 0 0.8])
cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Controllable_Sleep'
filename=strcat(num2str(time_shift),'_shift_Controllable_7am_1000lux','.mat');
load(filename);
hold on;
plot(x(:,1),H(:,2)-0.1333*x(:,2),'k','linewidth',6)
grid on
legend('Spontaneous sleep','Optimal sleep');
xlabel('time (hours)');ylabel('B');
box1=[-10 -10 300 300];
boy1=[0.17 0.27 0.27 0.17];
box2=[-10 -10 300 300];
boy2=[0.67 0.77 0.77 0.67];
hold on
patch(box1,boy1,[1 0 0],'FaceAlpha',0.2)
hold on
patch(box2,boy2,[1 0 0],'FaceAlpha',0.2)




% plot minimum-time entrainment strategy
clear all;
close all;
clc;
bdclose('all');
load('Periodic_Solution_JFK_I_1000lux.mat');
Imax=1000;
tol=0.01;Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
Initial_Time=1;
Periodic_Solution=Periodic_Solution(Initial_Time*100+1:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
load('12_shift_Controllable_7am_1000lux.mat');
I(:,2)=I(:,2).*(1-Sleep(:,2));
box1=[];box2=[];boy1=[];boy2=[];
if  I(1,2)>0
    box1 = [box1, 0   0];
    boy1 = [boy1, -1.5 1.5];
end
if  Sleep(1,2)>0
    box2 = [box2, 0   0];
    boy2 = [boy2, 0 1];
end
for i=2:size(x,1)
    if I(i,2)==0 && I(i-1,2)>0
        box1=[box1, I(i,1) I(i,1)];
        boy1=[boy1, 1.5 -1.5];
    elseif I(i,2)>0 && I(i-1,2)==0
        box1=[box1, I(i,1) I(i,1)];
        boy1=[boy1, -1.5 1.5];        
    end
    
    if Sleep(i,2)==0 && Sleep(i-1,2)>0
        box2=[box2, Sleep(i,1) Sleep(i,1)];
        boy2=[boy2, 1 0];
    elseif Sleep(i,2)>0 && Sleep(i-1,2)==0
        box2=[box2, Sleep(i,1) Sleep(i,1)];
        boy2=[boy2, 0 1];        
    end
    
    if i==size(x,1) && size(box1,2)/4~=round(size(box1,2)/4)
        if boy1(end)==1.5
            box1=[box1, I(i,1) I(i,1)];
            boy1=[boy1, 1.5 -1.5];
        elseif boy1(end)==0
            box1=[box1, I(i,1) I(i,1)];
            boy1=[boy1, -1.5 1.5];
        end
    end
    
    if i==size(x,1) && size(box2,2)/4~=round(size(box2,2)/4)
        if boy2(end)==1
            box2=[box2, Sleep(i,1) Sleep(i,1)];
            boy2=[boy2, 1 0];
        elseif boy2(end)==0
            box2=[box2, Sleep(i,1) Sleep(i,1)];
            boy2=[boy2, 0 1];
        end
    end   
end
subplot(2,1,1)
plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'b','linewidth',4)
hold on
plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',4)
hold on
plot(x(:,1),x(:,2),'k','linewidth',4)
hold on
plot(xc(:,1),xc(:,2),'k--','linewidth',4)
grid on
patch(box1,boy1,[1 0 0],'FaceAlpha',0.2)
legend('x_{ref}','x_{cref}','x','x_c','I*')
axis([0 x(end,1) -1.5 1.5])
subplot(2,1,2)
patch(box2,boy2,[0 0 0],'FaceAlpha',0.8)
hold on
plot(x(:,1),H(:,2)-0.1333*x(:,2),'b','linewidth',6)
grid on
legend('\beta','B')
hold on
axis([0 x(end,1) 0 1])
boX1=[-10 -10 300 300];
boY1=[0.17 0.27 0.27 0.17];
boX2=[-10 -10 300 300];
boY2=[0.67 0.77 0.77 0.67];
hold on
patch(boX1,boY1,[1 0 0],'FaceAlpha',0.2)
hold on
patch(boX2,boY2,[1 0 0],'FaceAlpha',0.2)