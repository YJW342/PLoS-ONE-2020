% Gradient descent method light of circadian entrainment of international travelers with spontaneous sleep;
clear all;
close all;
clc;
bdclose('all');
global Sleep u H x xc T_opt n I
load('Periodic_Solution_JFK_I_1000lux.mat');
load('f_Kronauer.mat');
Imax=1000;
tol=0.01;Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
Initial_Time=1;
Periodic_Solution=Periodic_Solution(Initial_Time*100+1:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
time_shift=[1:23];
for nn=1:size(time_shift,2)
    T_optimal=[];
    n0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),time_shift(nn));
    x0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),time_shift(nn));
    xc0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),time_shift(nn));
    H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),time_shift(nn));
    sim('Greedy_Advance_I.slx');
    T_Advance(nn)=x(end,1);
    sim('Greedy_Delay_I.slx');
    T_Delay(nn)=x(end,1);
    if T_Advance(nn)<T_Delay(nn)
        sim('Greedy_Advance_I.slx');
    end
    T_optimal(1)=x(end,1);T_opt=x(end,1)
    for j=1:50
        x_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),x(end,1));
        xc_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),x(end,1));
        H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),x(end,1));
        I_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),x(end,1));
        sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,6),x(end,1));
        n_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),x(end,1));
        u_ref_final=G*(1-n_ref_final)*alpha0*(I_ref_final/I0)^p*(1-sleep_ref_final);
        B_ref=(1-0.4*x_ref_final)*(1-0.4*xc_ref_final)*u_ref_final;
        dx_ref_final = (pi/12)*(xc_ref_final + B_ref  + mu*(1/3*x_ref_final + 4/3*x_ref_final^3 - 256/105*x_ref_final^7));
        dxc_ref_final = (pi/12)*(q*B_ref *xc_ref_final - (24/(taux*0.99729))^2*x_ref_final - k*B_ref*x_ref_final);
        dH_ref_final = (1-sleep_ref_final)*(1-H_ref_final)/tr-(H_ref_final*sleep_ref_final)/td;
        dn_ref_final = 60*(alpha0*(I_ref_final/I0)^p*(1-sleep_ref_final)*(1-n_ref_final)-beta*n_ref_final);
        Theta=(x(end,2)-x_ref_final)*(dx(end,2)-dx_ref_final)+(xc(end,2)-xc_ref_final)*(dxc(end,2)-dxc_ref_final)+(H(end,2)-H_ref_final)*(dH(end,2)-dH_ref_final);
        R1_final=-(x(end,2)-x_ref_final)/Theta;
        R2_final=-(xc(end,2)-xc_ref_final)/Theta;
        R3_final=-(H(end,2)-H_ref_final)/Theta;
        R4_final=0;          
        Switching_Time=[];
        for N=1:size(x,1)-1
            if Sleep(N,2)==1 && Sleep(N+1,2)==0
                Switching_Time=[Switching_Time,N];
            elseif Sleep(N,2)==0 && Sleep(N+1,2)==1
                Switching_Time=[Switching_Time,N];
            end    
        end
        if isempty(Switching_Time)
            tspan=[T_opt:-0.01:0];
            [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
            R=[fliplr(t');fliplr(r')]';
        else
            for N=1:size(Switching_Time,2)+1
                if N==1
                    tspan=[T_opt:-0.01:x(Switching_Time(size(Switching_Time,2)-N+1),1)+0.01];
                    [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
                    R=[fliplr(t');fliplr(r')]';
                elseif N==size(Switching_Time,2)+1
                    dot=dH(Switching_Time(size(Switching_Time,2)-N+2),2)-Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2);
                    R1_final=(dH(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,2) + dot*R(1,2) + Ac*dH(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,4) + Ac*dxc(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,3) + Ac*dn(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,5) - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,2) - Ac*dxc(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,3) - Ac*dH(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,4) - Ac*dn(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,5))/(dH(Switching_Time(size(Switching_Time,2)-N+2),2) + dot - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2));
                    R2_final=R(1,3);
                    R3_final=(dx(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,2) - dxc(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,3) - dn(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,5) - dx(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,2) + dxc(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,3) + dH(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,4) + dn(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,5) + dot*R(1,4) - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,4))/(dH(Switching_Time(size(Switching_Time,2)-N+2),2) + dot - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2));
                    R4_final=R(1,5);
                    tspan=[x(Switching_Time(size(Switching_Time,2)-N+2),1):-0.01:0];
                    [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
                    R=[[fliplr(t');fliplr(r')]';R];
                else
                    dot=dH(Switching_Time(size(Switching_Time,2)-N+2),2)-Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2);   
                    R1_final=(dH(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,2) + dot*R(1,2) + Ac*dH(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,4) + Ac*dxc(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,3) + Ac*dn(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,5) - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,2) - Ac*dxc(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,3) - Ac*dH(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,4) - Ac*dn(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,5))/(dH(Switching_Time(size(Switching_Time,2)-N+2),2) + dot - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2));
                    R2_final=R(1,3);
                    R3_final=(dx(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,2) - dxc(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,3) - dn(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,5) - dx(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,2) + dxc(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,3) + dH(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,4) + dn(Switching_Time(size(Switching_Time,2)-N+2)+1,2)*R(1,5) + dot*R(1,4) - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2)*R(1,4))/(dH(Switching_Time(size(Switching_Time,2)-N+2),2) + dot - Ac*dx(Switching_Time(size(Switching_Time,2)-N+2),2));
                    R4_final=R(1,5);
                    tspan=[x(Switching_Time(size(Switching_Time,2)-N+2),1):-0.01:x(Switching_Time(size(Switching_Time,2)-N+1),1)+0.01];
                    [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
                    R=[[fliplr(t');fliplr(r')]';R];
                end
            end
        end
        
        Gradient=pi*G*R(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(1-Sleep(:,2))/12+pi*G*R(:,3).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(q*xc(:,2)-k*x(:,2)).*(1-Sleep(:,2))/12+60*R(:,5).*(1-n(:,2)).*(1-Sleep(:,2));
        eta0=[0,10.^[-10:0.2:5]];
%         eta0=[0:0.01:1];
        for h=1:size(eta0,2)
            I_update=I;
            I_update(:,2)=I_update(:,2)-Gradient*eta0(h);
            I_update(:,2)=min(max(I_update(:,2),0),Imax);
            sim('Integration_Forward_JFK_Line_Search.slx');
            T0(h)=x_plus(end,1);
            Error0(h)=Error_plus(end,2);
        end
        [M,N]=min(T0);
        Num=find(T0==T0(N));
        [MM,NN]=min(Error0(Num));
        eta=eta0(Num(NN));
        if eta==0
            break
        else
            I_update=I;
            I_update(:,2)=I_update(:,2)-Gradient*eta;
            I_update(:,2)=min(max(I_update(:,2),0),Imax);
            sim('Integration_Forward_JFK.slx');
            T_opt=x(end,1)
            T_optimal(j+1)=T_opt;            
        end
    end
    filename=strcat(num2str(nn),'_shift_Spontaneous_Sleep_7am_1000lux','.mat');
    save(filename,'Sleep','x','xc','I','n','u','H','R','T_optimal','Gradient')
    T_Spontaneous(nn)=T_optimal(end);
    figure (nn)
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
    plot([0:size(Gradient,1)-1]*0.01,Gradient/max(abs(Gradient)),'b','linewidth',2)
    grid on
    axis([0 u(end,1) min(Gradient)/max(abs(Gradient)) 1])
end


    
    