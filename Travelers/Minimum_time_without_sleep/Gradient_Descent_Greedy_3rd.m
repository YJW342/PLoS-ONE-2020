% Golden Section Search
clear all
close all
clc
bdclose('all')
load('Periodic_Solution_JFK_1000lux.mat');
Periodic_Solution=Periodic_Solution(101:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
load('f_Kronauer.mat');
% xr=Periodic_Solution(:,2);
% xcr=Periodic_Solution(:,3);
time_shift=[1:23];
Imax=1000;
tol=0.01;Tol=0.05;
k=0.55;mu=0.13;q=1/3;tau=24.2;
G=33.75;I0=9500;p=0.5;alpha0=0.05;
for j=[12,14]%1:size(time_shift,2)
    T_OPT=[];
    x0=Periodic_Solution(time_shift(j)*100+1,2);
    xc0=Periodic_Solution(time_shift(j)*100+1,3);
    n0=Periodic_Solution(time_shift(j)*100+1,4);
    sim('Greedy_Advance_L.slx');
    T_Advance=x(end,1);
    sim('Greedy_Delay_L.slx');
    T_Delay=x(end,1);
    if T_Advance<T_Delay
        sim('Greedy_Advance_L.slx');
    end
%     I_update=Periodic_Solution(:,[1,6]);T_opt=400;
    I_update=I;T_opt=x(end,1);
    sim('Integration_Forward_Forger_L.slx');
    T_opt=x(end,1)
    for i=1:50
        % transversality condition
        xr_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),x(end,1));
        xcr_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),x(end,1));
        ur_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),x(end,1));
        B=(1-0.4*xr_final)*(1-0.4*xcr_final)*ur_final;
        dxr=pi*(xcr_final+mu*(xr_final/3+4*xr_final^3/3-256*xr_final^7/105)+B)/12;
        dxcr=pi*(q*B*xcr_final-xr_final*((24/0.99729/tau)^2+k*B))/12;
        dx_final=dx(end,2);dxc_final=dxc(end,2);   
        R1_final = -(x(end,2)-xr_final)/((x(end,2)-xr_final)*(dx_final-dxr)+(xc(end,2)-xcr_final)*(dxc_final-dxcr));
        R2_final = -(xc(end,2)-xcr_final)/((x(end,2)-xr_final)*(dx_final-dxr)+(xc(end,2)-xcr_final)*(dxc_final-dxcr));
        R3_final=0;
        T_opt=x(end,1)
        T_OPT(i)=T_opt;
        sim('Integration_Backward_Forger_L.slx');
        r1=fliplr(R1');r1=r1';
        r2=fliplr(R2');r2=r2';
        r3=fliplr(R3');r3=r3';
        R1(:,2)=r1(:,2);
        R2(:,2)=r2(:,2);
        R3(:,2)=r3(:,2);
        Gradient=(1-n(:,2)).*(pi*G*R1(:,2).*(1-0.4*x(:,2)).*(1-0.4*xc(:,2))/12+pi*R2(:,2)*G.*(1-0.4*x(:,2)).*(1-0.4*xc(:,2)).*[q*xc(:,2)-k*x(:,2)]/12+60*R3(:,2));
        eta0=[0,10.^[-5:0.2:5]];
        for h=1:size(eta0,2)
            I_update=I;
            I_update(:,2)=I_update(:,2)-eta0(h)*Gradient;
            I_update(:,2)=min(I_update(:,2),Imax);
            I_update(:,2)=max(I_update(:,2),0);
            sim('Integration_Forward_Forger_Line_Search.slx');
            T(h)=I_plus(end,1);
            error(h)=Error_plus(end,2);
        end
        [M,N]=min(T);
        Num=find(T==T(N));
        [MM,NN]=min(error(Num));
        eta=eta0(Num(NN));
        I_update=I;
        I_update(:,2)=I_update(:,2)-eta*Gradient;
        I_update(:,2)=min(I_update(:,2),Imax);
        I_update(:,2)=max(I_update(:,2),0);
        sim('Integration_Forward_Forger_L.mdl');
        if eta==0
            break
        end
    end

    filename=strcat(num2str(time_shift(j)),'_shift_1000lux','.mat');
    save(filename,'u','x','xc','n','R1','R2','R3','I','T_OPT')
    T_optimal(j)=T_OPT(end);
    figure(j)
    subplot(2,1,1)
    plot(Periodic_Solution(:,1),Periodic_Solution(:,2),'b','linewidth',2)
    hold on
    plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',2)
    hold on
    plot(x(:,1),x(:,2),'k','linewidth',2)
    hold on
    plot(x(:,1),xc(:,2),'k--','linewidth',2)
    hold on
    plot(u(:,1),u(:,2),'r','linewidth',2)
    axis([0 x(end,1) -1.3 1.2])
    grid on
    set(gca,'fontsize',16)
    xlabel('time/hours','FontSize', 16);ylabel('x','FontSize', 16);
    legend('x_{ref}','x_{cref}','x','x_c','u^*')
    subplot(2,1,2)
    plot([0:size(Gradient,1)-1]*0.01,Gradient/max(abs(Gradient)))
    hold on
    plot(I(:,1),I(:,2)/Imax)
    grid on
    xlabel('time(hours)')
    legend('\partialJ/\partialI','I^*')    
end


time_shift=[1:23];
for i=1:size(time_shift,2)
    filename=strcat(num2str(time_shift(i)),'_shift_1000lux','.mat');
    load(filename);
    T_optimal(i)=T_OPT(end);
end
hold on
plot(T_optimal)
