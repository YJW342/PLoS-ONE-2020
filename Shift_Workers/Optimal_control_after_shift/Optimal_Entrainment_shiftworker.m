% Optimal entrainment after shift with constant light during shift
clear all;
close all;
clc;
bdclose('all');
global Sleep u H x xc T_opt n I
load('Periodic_Solution_JFK_I_1000lux.mat');
Periodic_Solution=Periodic_Solution(1401:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
shift=12; % 12 represents shift work ends at 8 am;
Ishift=1000;Imax=10000;
tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
tol=0.01;
x0=Periodic_Solution(1,2); xc0=Periodic_Solution(1,3);
H0=Periodic_Solution(1,4); n0=Periodic_Solution(1,7);
Sleep0=0;
I_update=Periodic_Solution(:,[1,5]);
for i=1:size(I_update,1)
    if I_update(i,1)<shift %&& I_update(i,1)>2
        I_update(i,2)=Ishift;
    end
end
sim('Open_Loop_shift.slx');
T_optimal=[];
T_optimal(1)=x(end,1);T_opt=x(end,1)
x(end,1)-shift
for j=1:40
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
    dn_ref_final = 60*(alpha0*(1-sleep_ref_final)*(1-n_ref_final)*(I_ref_final/I0)^p-beta*n_ref_final);
    Theta=(x(end,2)-x_ref_final)*(dx(end,2)-dx_ref_final)+(xc(end,2)-xc_ref_final)*(dxc(end,2)-dxc_ref_final)+(H(end,2)-H_ref_final)*(dH(end,2)-dH_ref_final);
    R1_final=-(x(end,2)-x_ref_final)/Theta;
    R2_final=-(xc(end,2)-xc_ref_final)/Theta;
    R3_final=-(H(end,2)-H_ref_final)/Theta;
    R4_final=0;
    tspan=[T_opt:-0.01:0];
    [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
    R=[fliplr(t');fliplr(r')]';     
    Gradient=pi*G*R(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(1-Sleep(:,2))/12+pi*G*R(:,3).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(q*xc(:,2)-k*x(:,2)).*(1-Sleep(:,2))/12+60*R(:,5).*(1-n(:,2)).*(1-Sleep(:,2));
    eta0=[0,10.^[-10:0.2:2]];
    Sleep_update=Sleep;
    for i=1:size(x,1)
        if x(i,1)<shift
            Gradient(i)=0;
        end
    end
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
    eta1=eta0(Num(NN));
    I_update=I;
    I_update(:,2)=I_update(:,2)-Gradient*eta1;
    I_update(:,2)=min(max(I_update(:,2),0),Imax);
    sim('Integration_Forward_JFK.slx');
    T_opt=x(end,1)
    T_optimal=[T_optimal,T_opt];
    x_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),x(end,1));
    xc_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),x(end,1));
    H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),x(end,1));
    I_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),x(end,1));
    sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,6),x(end,1));
    n_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),x(end,1));
    u_ref_final=G*(1-n_ref_final)*alpha0*(1-sleep_ref_final)*(I_ref_final/I0)^p;       
    B_ref=(1-0.4*x_ref_final)*(1-0.4*xc_ref_final)*u_ref_final;
    dx_ref_final = (pi/12)*(xc_ref_final + B_ref  + mu*(1/3*x_ref_final + 4/3*x_ref_final^3 - 256/105*x_ref_final^7));
    dxc_ref_final = (pi/12)*(q*B_ref *xc_ref_final - (24/(taux*0.99729))^2*x_ref_final - k*B_ref*x_ref_final);
    dH_ref_final = (1-sleep_ref_final)*(1-H_ref_final)/tr-(H_ref_final*sleep_ref_final)/td;
    dn_ref_final = 60*(alpha0*(1-n_ref_final)*(1-sleep_ref_final)*(I_ref_final/I0)^p-beta*n_ref_final);
    Theta=(x(end,2)-x_ref_final)*(dx(end,2)-dx_ref_final)+(xc(end,2)-xc_ref_final)*(dxc(end,2)-dxc_ref_final)+(H(end,2)-H_ref_final)*(dH(end,2)-dH_ref_final);
    R1_final=-(x(end,2)-x_ref_final)/Theta;
    R2_final=-(xc(end,2)-xc_ref_final)/Theta;
    R3_final=-(H(end,2)-H_ref_final)/Theta;
    R4_final=0;
    tspan=[T_opt:-0.01:0];
    [t,r]=ode113(@(t,R) costate(t,R),tspan, [R1_final;R2_final;R3_final;R4_final]);
    R=[fliplr(t');fliplr(r')]';
    Gradient=pi*G*R(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(1-Sleep(:,2))/12+pi*G*R(:,3).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*(1-n(:,2)).*(q*xc(:,2)-k*x(:,2)).*(1-Sleep(:,2))/12+60*R(:,5).*(1-n(:,2)).*(1-Sleep(:,2));
    
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
        delta_switching(i)=R(Switching_Time(i),2)*[dx(Switching_Time(i),2)-dx(Switching_Time(i)+1,2)]+R(Switching_Time(i),3)*[dxc(Switching_Time(i),2)-dxc(Switching_Time(i)+1,2)]+R(Switching_Time(i),4)*[dH(Switching_Time(i),2)-dH(Switching_Time(i)+1,2)]+R(Switching_Time(i),5)*[dn(Switching_Time(i),2)-dn(Switching_Time(i)+1,2)];
    end
        
    eta0=[0:50];
    I_update=I;
    for h=1:size(eta0,2)
        Sleep_update=Sleep;
        for i=2:size(Switching_Time,2)
            if delta_switching(i)>0
                Sleep_update(max([Switching_Time(i)+1-round(eta0(h)*delta_switching(i)/abs(delta_switching(i))),1]):Switching_Time(i)+1,2)=Sleep_update(Switching_Time(i)+1,2);
            elseif delta_switching(i)<0
                Sleep_update(Switching_Time(i):Switching_Time(i)-round(eta0(h)*delta_switching(i)/abs(delta_switching(i))),2)=Sleep_update(Switching_Time(i),2);
            end
        end
        sim('Integration_Forward_JFK_Line_Search.slx')
        T(h)=x_plus(end,1);
        error(h)=Error_plus(end,2);
    end
    [M,N]=min(T);
    Num=find(T==T(N));
    [MM,NN]=min(error(Num));
    eta2=eta0(Num(NN));
    Sleep_update=Sleep;
    for i=2:size(Switching_Time,2)
        if delta_switching(i)>0
            Sleep_update(max([Switching_Time(i)+1-round(eta2*delta_switching(i)/abs(delta_switching(i))),1]):Switching_Time(i)+1,2)=Sleep_update(Switching_Time(i)+1,2);
        elseif delta_switching(i)<0
            Sleep_update(Switching_Time(i):Switching_Time(i)-round(eta2*delta_switching(i)/abs(delta_switching(i))),2)=Sleep_update(Switching_Time(i),2);
        end
    end
    sim('Integration_Forward_JFK.slx')
    T_opt=x(end,1)
    T_optimal=[T_optimal,T_opt];
    
    if eta1==0 && eta2==0
        break
    end
end
filename=strcat('shift_worker_0lux_10000lux_10pm','.mat');
save(filename,'Sleep','x','xc','I','u','H','R','n','T_optimal')
