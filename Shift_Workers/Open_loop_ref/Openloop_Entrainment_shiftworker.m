% Optimal entrainment during shift with open-loop entrainment afterware
clear all;
% close all;
clc;
bdclose('all');
load('Periodic_Solution_JFK_I_1000lux.mat');
Periodic_Solution=Periodic_Solution(1601:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
shift=10; % 12 represents shift work ends at 8 am;
Ishift=1000;
tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
tol=0.01;
x0=Periodic_Solution(1,2); xc0=Periodic_Solution(1,3);
H0=Periodic_Solution(1,4); n0=Periodic_Solution(1,7);
Sleep0=0;
% I_update=Periodic_Solution(:,[1,5]);
% for i=1:size(I_update,1)
%     if I_update(i,1)<shift %&& I_update(i,1)>2
%         I_update(i,2)=Ishift;
%     end
% end
sim('Open_Loop_shift.slx');
x(end,1)-shift


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plot schedule
clear all;
% close all;
clc;
bdclose('all');
load('Periodic_Solution_JFK_I_1000lux.mat');
Periodic_Solution=Periodic_Solution(1601:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
shift=10; % 12 represents shift work ends at 8 am;
Ishift=0;
tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
tol=0.01;
x0=Periodic_Solution(1,2); xc0=Periodic_Solution(1,3);
H0=Periodic_Solution(1,4); n0=Periodic_Solution(1,7);
I_update=Periodic_Solution(:,[1,5]);
for i=1:size(I_update,1)
    if I_update(i,1)<shift %&& I_update(i,1)>2
        I_update(i,2)=Ishift;
    end
end
sim('Open_Loop_shift.slx');
I(:,2)=I(:,2).*(1-Sleep(:,2));
box1=[];box2=[];boy1=[];boy2=[];
if  I(1,2)>0
    box1 = [box1, 0   0];
    boy1 = [boy1, 0 1.5];
end
if  Periodic_Solution(1,5)>0
    box2 = [box2, 0   0];
    boy2 = [boy2, 0 -1.5];
end

for i=2:size(x,1)
    if I(i,2)==0 && I(i-1,2)>0
        box1=[box1, I(i,1) I(i,1)];
        boy1=[boy1, 1.5 0];
    elseif I(i,2)>0 && I(i-1,2)==0
        box1=[box1, I(i,1) I(i,1)];
        boy1=[boy1, 0 1.5];        
    end
    
    if Periodic_Solution(i,5)==0 && Periodic_Solution(i-1,5)>0
        box2=[box2, Periodic_Solution(i,1) Periodic_Solution(i,1)];
        boy2=[boy2, -1.5 0];
    elseif Periodic_Solution(i,5)>0 && Periodic_Solution(i-1,5)==0
        box2=[box2, Periodic_Solution(i,1) Periodic_Solution(i,1)];
        boy2=[boy2, 0 -1.5];        
    end
    if i==size(x,1) && size(box1,2)/4~=round(size(box1,2)/4)
        if boy1(end)==1.5
            box1=[box1, I(i,1) I(i,1)];
            boy1=[boy1, 1.5 0];
        elseif boy1(end)==0
            box1=[box1, I(i,1) I(i,1)];
            boy1=[boy1, 0 1.5];
        end
    end
    if i==size(x,1) && size(box2,2)/4~=round(size(box2,2)/4)
        if boy2(end)==-1.5
            box2=[box2, I(i,1) I(i,1)];
            boy2=[boy2, -1.5 0];
        elseif boy2(end)==0
            box2=[box2, I(i,1) I(i,1)];
            boy2=[boy2, 0 -1.5];
        end
    end
end


box3=[];box4=[];boy3=[];boy4=[];
if  Sleep(1,2)>0
    box3 = [box3, 0   0];
    boy3 = [boy3, 0 0.5];
end
if  Periodic_Solution(1,6)>0
    box4 = [box4, 0   0];
    boy4 = [boy4, 0 -0.5];
end

for i=2:size(x,1)
    if Sleep(i,2)==0 && Sleep(i-1,2)>0
        box3=[box3, Sleep(i,1) Sleep(i,1)];
        boy3=[boy3, .5 0];
    elseif Sleep(i,2)>0 && Sleep(i-1,2)==0
        box3=[box3, Sleep(i,1) Sleep(i,1)];
        boy3=[boy3, 0 .5];        
    end
    
    if Periodic_Solution(i,6)==0 && Periodic_Solution(i-1,6)>0
        box4=[box4, Periodic_Solution(i,1) Periodic_Solution(i,1)];
        boy4=[boy4, -.5 0];
    elseif Periodic_Solution(i,6)>0 && Periodic_Solution(i-1,6)==0
        box4=[box4, Periodic_Solution(i,1) Periodic_Solution(i,1)];
        boy4=[boy4, 0 -.5];        
    end
    if i==size(x,1) && size(box3,2)/4~=round(size(box3,2)/4)
        if boy3(end)==0.5
            box3=[box3, Sleep(i,1) Sleep(i,1)];
            boy3=[boy3, 0.5 0];
        elseif boy3(end)==0
            box3=[box3, Sleep(i,1) Sleep(i,1)];
            boy3=[boy3, 0 0.5];
        end
    end
    if i==size(x,1) && size(box4,2)/4~=round(size(box4,2)/4)
        if boy4(end)==-0.5
            box4=[box4, Sleep(i,1) Sleep(i,1)];
            boy4=[boy4, -0.5 0];
        elseif boy4(end)==0
            box4=[box4, Sleep(i,1) Sleep(i,1)];
            boy4=[boy4, 0 -0.5];
        end
    end
end


subplot(3,1,3)
patch(box1-shift,boy1,[1 0 0],'FaceAlpha',0.8)
patch(box2-shift,boy2,[1 0 0],'FaceAlpha',0.2)
patch(box3-shift,2*boy3,[0 0 0],'FaceAlpha',0.8)
patch(box4-shift,2*boy4,[0 0 0],'FaceAlpha',0.2)
hold on
plot(x(:,1)-shift,x(:,2),'b','linewidth',6)
hold on
plot(xc(:,1)-shift,xc(:,2),'g','linewidth',6)
hold on
plot(H(:,1)-shift,H(:,2),'y','linewidth',6)
hold on
% plot(Sleep(:,1),Sleep(:,2),'k','linewidth',3)
% hold on
plot(Periodic_Solution(:,1)-shift,Periodic_Solution(:,2),'b--','linewidth',6)
hold on
plot(Periodic_Solution(:,1)-shift,Periodic_Solution(:,3),'g--','linewidth',6)
hold on
plot(Periodic_Solution(:,1)-shift,Periodic_Solution(:,4),'y--','linewidth',6)
hold on
% plot(Periodic_Solution(:,1),Periodic_Solution(:,6),'k--','linewidth',3)
% xlabel('time (h)')
grid on
axis([-shift x(end,1)-shift -1.3 1.2])
legend('I^*','I_{ref}','\beta^*','\beta_{ref}','x','x_c','H','x_{ref}','x_{cref}','H_{ref}')
xlabel('time (h)')
