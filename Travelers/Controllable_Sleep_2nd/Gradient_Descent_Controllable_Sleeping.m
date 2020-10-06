%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Jiawei Yin
% 02.27.2020
% Gradient descent method for sleep schedule and light of circadian entrainment of international travelers on S+C2 model;
clear all;
close all;
clc;
bdclose('all');
global Sleep u H x xc T_opt
load('Periodic_Solution_JFK_u_1000lux.mat');
load('f_Kronauer.mat');
umax=0.2208;
tol=0.01;Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;
Initial_Time=1;
Periodic_Solution=Periodic_Solution(Initial_Time*100+1:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
time_shift=[1:23];
for n=1:size(time_shift,2)
    T_optimal=[];
    x0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),time_shift(n));
    xc0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),time_shift(n));
    H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),time_shift(n));
    Sleep0=0;
    % give initial guess of u and sleeping and waking time
    sim('Greedy_Advance_u.slx');
    T_Advance(n)=x(end,1);
    sim('Greedy_Delay_u.slx');
    T_Delay(n)=x(end,1);
    if T_Advance(n)<T_Delay(n)
        sim('Greedy_Advance_u.slx');
    end
    T_optimal(1)=x(end,1);T_opt=x(end,1);
    for j=1:30
        x_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),x(end,1));
        xc_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),x(end,1));
        H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),x(end,1));
        u_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),x(end,1));
        sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,6),x(end,1));
        B_ref=(1-0.4*x_ref_final)*(1-0.4*xc_ref_final)*u_ref_final*(1-sleep_ref_final);
        dx_ref_final = (pi/12)*(xc_ref_final + B_ref  + mu*(1/3*x_ref_final + 4/3*x_ref_final^3 - 256/105*x_ref_final^7));
        dxc_ref_final = (pi/12)*(q*B_ref *xc_ref_final - (24/(taux*0.99729))^2*x_ref_final - k*B_ref*x_ref_final);
        dH_ref_final = (1-sleep_ref_final)*(1-H_ref_final)/tr-(H_ref_final*sleep_ref_final)/td;
        Theta=(x(end,2)-x_ref_final)*(dx(end,2)-dx_ref_final)+(xc(end,2)-xc_ref_final)*(dxc(end,2)-dxc_ref_final)+(H(end,2)-H_ref_final)*(dH(end,2)-dH_ref_final);
        R1_final=-(x(end,2)-x_ref_final)/Theta;
        R2_final=-(xc(end,2)-xc_ref_final)/Theta;
        R3_final=-(H(end,2)-H_ref_final)/Theta;
        tspan=[T_opt:-0.01:0];
        [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final]);
        R=[fliplr(t');fliplr(r')]';     
        Gradient=pi*R(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-Sleep(:,2))/12+pi*R(:,3).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-Sleep(:,2)).*(q*xc(:,2)-k*x(:,2))/12;
%         eta0=[0,10^[-5:0.1:0]];
        eta0=[0,10.^[-10:0.2:0]];
        Sleep_update=Sleep;
        for h=1:size(eta0,2)
            u_update=u;
            u_update(:,2)=min(max(u(:,2)-eta0(h)*Gradient,0),umax);
            sim('Integration_Forward_JFK_Line_Search.slx');
            T0(h)=x_plus(end,1);
            Error0(h)=Error_plus(end,2);
        end
        [M,N]=min(T0);
        Num=find(T0==T0(N));
        [MM,NN]=min(Error0(Num));
        eta1=eta0(Num(NN));
        u_update=u;
        u_update(:,2)=min(max(u(:,2)-eta1*Gradient,0),umax);
        sim('Integration_Forward_JFK.slx');
%         T_opt=x(end,1);
        x_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),x(end,1));
        xc_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),x(end,1));
        H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),x(end,1));
        u_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),x(end,1));
        sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,6),x(end,1));
        B_ref=(1-0.4*x_ref_final)*(1-0.4*xc_ref_final)*u_ref_final*(1-sleep_ref_final);
        dx_ref_final = (pi/12)*(xc_ref_final + B_ref  + mu*(1/3*x_ref_final + 4/3*x_ref_final^3 - 256/105*x_ref_final^7));
        dxc_ref_final = (pi/12)*(q*B_ref *xc_ref_final - (24/(taux*0.99729))^2*x_ref_final - k*B_ref*x_ref_final);
        dH_ref_final = (1-sleep_ref_final)*(1-H_ref_final)/tr-(H_ref_final*sleep_ref_final)/td;
        Theta=(x(end,2)-x_ref_final)*(dx(end,2)-dx_ref_final)+(xc(end,2)-xc_ref_final)*(dxc(end,2)-dxc_ref_final)+(H(end,2)-H_ref_final)*(dH(end,2)-dH_ref_final);
        R1_final=-(x(end,2)-x_ref_final)/Theta;
        R2_final=-(xc(end,2)-xc_ref_final)/Theta;
        R3_final=-(H(end,2)-H_ref_final)/Theta;
        T_opt=x(end,1)
        tspan=[T_opt:-0.01:0];
        [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final]);
        R=[fliplr(t');fliplr(r')]';
        Gradient=pi*R(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-Sleep(:,2))/12+pi*R(:,3).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-Sleep(:,2)).*(q*xc(:,2)-k*x(:,2))/12;
        
        Switching_Time=[];
        for N=1:size(x,1)-1
            if Sleep(N,2)==1 && Sleep(N+1,2)==0
                Switching_Time=[Switching_Time,N];
            elseif Sleep(N,2)==0 && Sleep(N+1,2)==1
                Switching_Time=[Switching_Time,N];
            end    
        end
        
        delta_switching=[];
        for i=1:size(Switching_Time,2)
            delta_switching(i)=R(Switching_Time(i),2)*[dx(Switching_Time(i),2)-dx(Switching_Time(i)+1,2)]+R(Switching_Time(i),3)*[dxc(Switching_Time(i),2)-dxc(Switching_Time(i)+1,2)]+R(Switching_Time(i),4)*[dH(Switching_Time(i),2)-dH(Switching_Time(i)+1,2)];
        end
        
        eta0=[0:50];
        u_update=u;
        for h=1:size(eta0,2)
            Sleep_update=Sleep;
            for i=1:size(Switching_Time,2)
                if delta_switching(i)>0
                    Sleep_update(max([Switching_Time(i)+1-round(eta0(h)*delta_switching(i)/abs(delta_switching(i))),1]):Switching_Time(i)+1,2)=Sleep_update(Switching_Time(i)+1,2);
                elseif delta_switching(i)<0
                    Sleep_update(Switching_Time(i):Switching_Time(i)-round(eta0(h)*delta_switching(i)/abs(delta_switching(i))),2)=Sleep_update(Switching_Time(i),2);
                end
            end
            sim('Gradient_Descent_Sleep_Line_Search.slx')
            T(h)=x_plus(end,1);
            error(h)=Error_plus(end,2);
        end
        [M,N]=min(T);
        Num=find(T==T(N));
        [MM,NN]=min(error(Num));
        eta2=eta0(Num(NN));
        Sleep_update=Sleep;
        for i=1:size(Switching_Time,2)
            if delta_switching(i)>0
                Sleep_update(max([Switching_Time(i)+1-round(eta2*delta_switching(i)/abs(delta_switching(i))),1]):Switching_Time(i)+1,2)=Sleep_update(Switching_Time(i)+1,2);
            elseif delta_switching(i)<0
                Sleep_update(Switching_Time(i):Switching_Time(i)-round(eta2*delta_switching(i)/abs(delta_switching(i))),2)=Sleep_update(Switching_Time(i),2);
            end
        end
        if eta1==0 && eta2==0
            break
        else
            sim('Gradient_Descent_Sleep.slx')
            T_opt=x(end,1)
            T_optimal(j+1)=T_opt;
        end
    end
    filename=strcat(num2str(n),'_shift_Controllable_Sleep_7am_10000lux','.mat');
    save(filename,'Sleep','x','xc','u','H','R','T_optimal','Gradient')
    T_Controllable(n)=T_optimal(end);
    figure (n)
    subplot(2,1,1)
    plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'b','linewidth',2)
    hold on
    plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',2)
    hold on
    plot(x(:,1),x(:,2),'k','linewidth',2)
    hold on
    plot(xc(:,1),xc(:,2),'k--','linewidth',2)
    hold on
    plot(u(:,1),u(:,2),'r','linewidth',2)
    grid on
    axis([0 u(end,1) -1.3 1.3])
    subplot(2,1,2)
    plot(u(:,1),u(:,2),'r','linewidth',2)
    hold on
    plot([0:size(Gradient,1)-1]*0.01,Gradient/100,'b','linewidth',2)
    hold on
    plot(B(:,1),B(:,2),'k','linewidth',2)
    grid on
end
figure (24)
plot(T_Controllable)
grid on