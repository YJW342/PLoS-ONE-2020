% Gradient descent method for sleep schedule and light of circadian entrainment of international travelers;
% fixed reference condition
clear all;
close all;
clc;
bdclose('all');
global Sleep u H f theta T_opt
load('Periodic_Solution_PRC_u_1000lux.mat');
Periodic_Solution=Periodic_Solution(101:end,:);Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
load('f_Kronauer.mat');
umax=0.2208;
tol=0.01;time_shift=[1:23];
for n=1:size(time_shift,2)
    T_optimal=[];
    theta0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),time_shift(n));
    H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),time_shift(n));
    Sleep0=0;
    sim('Greedy_Advance_u.slx');
    T_Advance(n)=theta(end,1);
    sim('Greedy_Delay_u.slx');
    T_Delay(n)=theta(end,1);
    if T_Advance(n)<T_Delay(n)
        sim('Greedy_Advance_u.slx');
    end
    T_Greedy(n)=theta(end,1);
    T_opt=theta(end,1)
    T_optimal(1)=theta(end,1);T_opt=theta(end,1);
    for j=1:30
        theta_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),theta(end,1));
        H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),theta(end,1));
        u_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),theta(end,1));
        sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),theta(end,1));
        theta_diff=[(theta(end,2)-theta_ref_final),(theta(end,2)-theta_ref_final+2*pi),(theta(end,2)-theta_ref_final-2*pi)];
        [M,N]=min(abs(theta_diff));
        Theta_diff=theta_diff(N);
        R1_final=-Theta_diff/[Theta_diff*(interp1q(f(:,1),f(:,2),theta(end,2))*u(end,2)*(1-Sleep(end,2))-interp1q(f(:,1),f(:,2),theta_ref_final)*(1-sleep_ref_final)*u_ref_final)+(H(end,2)-H_ref_final)*((1-Sleep(end,2))*(1-H(end,2))/18.2-Sleep(end,2)*H(end,2)/4.2-(1-sleep_ref_final)*(1-H_ref_final)/18.2+sleep_ref_final*H_ref_final/4.2)];
        R2_final=-(H(end,2)-H_ref_final)/[Theta_diff*(interp1q(f(:,1),f(:,2),theta(end,2))*u(end,2)*(1-Sleep(end,2))-interp1q(f(:,1),f(:,2),theta_ref_final)*(1-sleep_ref_final)*u_ref_final)+(H(end,2)-H_ref_final)*((1-Sleep(end,2))*(1-H(end,2))/18.2-Sleep(end,2)*H(end,2)/4.2-(1-sleep_ref_final)*(1-H_ref_final)/18.2+sleep_ref_final*H_ref_final/4.2)];
        tspan=[T_opt:-0.01:0];
        [t,r1]=ode113(@(t,R1) costate1(t,R1),tspan, R1_final);
        R1=[fliplr(t');fliplr(r1')]';
        [t,r2]=ode113(@(t,R2) costate2(t,R2),tspan, R2_final);
        R2=[fliplr(t');fliplr(r2')]';
        Waking_Time=[];Sleeping_Time=[];
        for N=1:size(theta,1)-1
            if Sleep(N,2)==1 && Sleep(N+1,2)==0
                Waking_Time=[Waking_Time,N];
            end
            if Sleep(N,2)==0 && Sleep(N+1,2)==1
                Sleeping_Time=[Sleeping_Time,N];
            end
        end
        delta_waking=[];delta_sleeping=[];
        for i=1:size(Waking_Time,2)
            delta_waking(i)=-R1(Waking_Time(i),2)*interp1q(f(:,1),f(:,2),theta(Waking_Time(i)+1,2))*u(Waking_Time(i)+1,2)-R2(Waking_Time(i),2)*((1-H(Waking_Time(i)+1,2))/18.2+H(Waking_Time(i),2)/4.2);
        end
        for i=1:size(Sleeping_Time,2)
            delta_sleeping(i)=R1(Sleeping_Time(i),2)*interp1q(f(:,1),f(:,2),theta(Sleeping_Time(i),2))*u(Sleeping_Time(i),2)+R2(Sleeping_Time(i),2)*((1-H(Sleeping_Time(i),2))/18.2+H(Sleeping_Time(i)+1,2)/4.2);
        end
        % delta_waking
        % delta_sleeping
        % line_search
        eta0=[0:50];
        if R1(end,2)<0
            for h=1:size(eta0,2)
                Sleep_update=Sleep;
                for i=1:size(Waking_Time,2)
                    if delta_waking(i)>0
                        Sleep_update(max([Waking_Time(i)+1-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),1]):Waking_Time(i)+1,2)=0;
                    elseif delta_waking(i)<0
                        Sleep_update(Waking_Time(i):Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),2)=1;
                    end
                end
                for i=1:size(Sleeping_Time,2)
                    if delta_sleeping(i)<0
                        Sleep_update(Sleeping_Time(i):Sleeping_Time(i)-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                    elseif delta_sleeping(i)>0
                        Sleep_update(max([Sleeping_Time(i)+1-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),1]):Sleeping_Time(i)+1,2)=1;
                    end
                end
                sim('Gradient_Descent_Sleep_Advance_Line_Search.slx')
                T(h)=theta_plus(end,1);
                error(h)=Error(end,2);
            end
            [M,N]=min(T);
            Num=find(T==T(N));
            [MM,NN]=min(error(Num));
            eta=eta0(Num(NN));
            for i=1:size(Waking_Time,2)
                if delta_waking(i)>0
                    Sleep(max([Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))),1]):Waking_Time(i),2)=0;
                elseif delta_waking(i)<0
                    Sleep(Waking_Time(i):Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))),2)=1;
                end
            end
            for i=1:size(Sleeping_Time,2)
                if delta_sleeping(i)<0
                    Sleep(Sleeping_Time(i):Sleeping_Time(i)-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                elseif delta_sleeping(i)>0
                    Sleep(max([Sleeping_Time(i)+1-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))),1]):Sleeping_Time(i)+1,2)=1;
                end
            end
            sim('Gradient_Descent_Sleep_Advance.slx')
            T_opt=theta(end,1)
        else
            for h=1:size(eta0,2)
                Sleep_update=Sleep;
                for i=1:size(Waking_Time,2)
                    if delta_waking(i)>0
                        Sleep_update(max([Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),1]):Waking_Time(i),2)=0;
                    elseif delta_waking(i)<0
                        Sleep_update(Waking_Time(i):Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),2)=1;
                    end
                end
                for i=1:size(Sleeping_Time,2)
                    if delta_sleeping(i)<0
                        Sleep_update(Sleeping_Time(i):Sleeping_Time(i)-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                    elseif delta_sleeping(i)>0
                        Sleep_update(max([Sleeping_Time(i)+1-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),1]):Sleeping_Time(i)+1,2)=1;
                    end
                end
                sim('Gradient_Descent_Sleep_Delay_Line_Search.slx')
                T(h)=theta_plus(end,1);
                error(h)=Error(end,2);
            end
            [M,N]=min(T);
            Num=find(T==T(N));
            [MM,NN]=min(error(Num));
            eta=eta0(Num(NN));
            for i=1:size(Waking_Time,2)
                if delta_waking(i)>0
                    Sleep(Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))):Waking_Time(i),2)=0;
                elseif delta_waking(i)<0
                    Sleep(Waking_Time(i):Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))),2)=1;
                end
            end
            for i=1:size(Sleeping_Time,2)
                if delta_sleeping(i)<0
                    Sleep(Sleeping_Time(i):Sleeping_Time(i)-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                elseif delta_sleeping(i)>0
                    Sleep(Sleeping_Time(i)+1-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))):Sleeping_Time(i)+1,2)=1;
                end
            end
            sim('Gradient_Descent_Sleep_Delay.slx')
            T_opt=theta(end,1)
        end
        %  clear theta u H
        T_optimal(j+1)=T_opt;
        Sleep=sleep;
        if eta==0
            break
        end
    end
    N=min([size(R1,1),size(theta,1)]);
    filename=strcat(num2str(n),'_shift_Controllable_Sleep_7am_10000lux','.mat');
    save(filename,'sleep','theta','u','H','R1','R2','T_optimal')
    T_Controllable(n)=T_optimal(end);
    figure (n)
    subplot(2,1,1)
    plot(Periodic_Solution(:,1),Periodic_Solution(:,2)/2/pi,'b','linewidth',2)
    hold on
    plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',2)
    hold on
    plot(theta(:,1),theta(:,2)/2/pi,'k','linewidth',2)
    hold on
    plot(H(:,1),H(:,2),'k--','linewidth',2)
    hold on
    plot(u(:,1),u(:,2),'r','linewidth',2)
    grid on
    xlabel('time(hours)');
    axis([0 theta(end,1) 0 1.2])
    legend('\theta_{ref}','H_{ref}','\theta','H','u')
    set(gca,'fontsize',20)
    subplot(2,1,2)
    plot(R1(1:N,1),R1(1:N,2).*(1-sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20,'b','linewidth',2)
    hold on
    plot(u(:,1),u(:,2).*(1-sleep(:,2)),'r','linewidth',2)
    hold on
    plot(Sleep(:,1),Sleep(:,2),'k','linewidth',2)
    grid on
    xlabel('time(hours)');
    axis([0 theta(end,1) min(R1(1:N,2).*(1-sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20) max(R1(1:N,2).*(1-Sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20)])
%     legend('P_1(1-\beta)f(\theta)','u')
    set(gca,'fontsize',20)
end




% Gradient descent method for sleep schedule and light of circadian entrainment of international travelers;
% fixed initial condition
clear all;
clc;
bdclose('all');
global Sleep u H f theta T_opt
load('Periodic_Solution_PRC_sleep.mat');
load('f_Kronauer.mat');
umax=0.2208;
tol=0.01;
Initial_Time=21;% initial local time of subjects
theta0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),Initial_Time-6);
H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),Initial_Time-6);
time_shift=[1:23];
Sleep0=0;
for n=1:size(time_shift,2)
    T_optimal=[];
    load('Periodic_Solution_PRC_sleep.mat');
    Time_shift=time_shift(n);% time shift for international travelers
    Initial_Reference=mod(Initial_Time+Time_shift-6,24);
    [M,N]=min(abs([Periodic_Solution(:,1)-Initial_Reference]))
    Periodic_Solution=Periodic_Solution(N:end,:);
    Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
    % give initial guess of u and sleeping and waking time
    sim('Greedy_Advance_u.slx');
    T_Advance(n)=theta(end,1);
    sim('Greedy_Delay_u.slx');
    T_Delay(n)=theta(end,1);
    if T_Advance(n)<T_Delay(n)
        sim('Greedy_Advance_u.slx');
    end
    T_Greedy(n)=theta(end,1);
    T_opt=theta(end,1)
    T_optimal(1)=theta(end,1);T_opt=theta(end,1);
    for j=1:50
        theta_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),theta(end,1));
        H_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),theta(end,1));
        u_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),theta(end,1));
        sleep_ref_final=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,5),theta(end,1));
        theta_diff=[(theta(end,2)-theta_ref_final),(theta(end,2)-theta_ref_final+2*pi),(theta(end,2)-theta_ref_final-2*pi)];
        [M,N]=min(abs(theta_diff));
        Theta_diff=theta_diff(N);
        R1_final=-Theta_diff/[Theta_diff*(interp1q(f(:,1),f(:,2),theta(end,2))*u(end,2)*(1-Sleep(end,2))-interp1q(f(:,1),f(:,2),theta_ref_final)*(1-sleep_ref_final)*u_ref_final)+(H(end,2)-H_ref_final)*((1-Sleep(end,2))*(1-H(end,2))/18.2-Sleep(end,2)*H(end,2)/4.2-(1-sleep_ref_final)*(1-H_ref_final)/18.2+sleep_ref_final*H_ref_final/4.2)];
        R2_final=-(H(end,2)-H_ref_final)/[Theta_diff*(interp1q(f(:,1),f(:,2),theta(end,2))*u(end,2)*(1-Sleep(end,2))-interp1q(f(:,1),f(:,2),theta_ref_final)*(1-sleep_ref_final)*u_ref_final)+(H(end,2)-H_ref_final)*((1-Sleep(end,2))*(1-H(end,2))/18.2-Sleep(end,2)*H(end,2)/4.2-(1-sleep_ref_final)*(1-H_ref_final)/18.2+sleep_ref_final*H_ref_final/4.2)];
        tspan=[T_opt:-0.01:0];
        [t,r1]=ode113(@(t,R1) costate1(t,R1),tspan, R1_final);
        R1=[fliplr(t');fliplr(r1')]';
        [t,r2]=ode113(@(t,R2) costate2(t,R2),tspan, R2_final);
        R2=[fliplr(t');fliplr(r2')]';
        Waking_Time=[];Sleeping_Time=[];
        for N=1:size(theta,1)-1
            if Sleep(N,2)==1 && Sleep(N+1,2)==0
                Waking_Time=[Waking_Time,N];
            end
            if Sleep(N,2)==0 && Sleep(N+1,2)==1
                Sleeping_Time=[Sleeping_Time,N];
            end
        end
        delta_waking=[];delta_sleeping=[];
        for i=1:size(Waking_Time,2)
            delta_waking(i)=-R1(Waking_Time(i),2)*interp1q(f(:,1),f(:,2),theta(Waking_Time(i)+1,2))*u(Waking_Time(i)+1,2)-R2(Waking_Time(i),2)*((1-H(Waking_Time(i)+1,2))/18.2+H(Waking_Time(i),2)/4.2);
        end
        for i=1:size(Sleeping_Time,2)
            delta_sleeping(i)=R1(Sleeping_Time(i),2)*interp1q(f(:,1),f(:,2),theta(Sleeping_Time(i),2))*u(Sleeping_Time(i),2)+R2(Sleeping_Time(i),2)*((1-H(Sleeping_Time(i),2))/18.2+H(Sleeping_Time(i)+1,2)/4.2);
        end
        % delta_waking
        % delta_sleeping
        % line_search
        eta0=[0:100];
        if R1(end,2)<0
            for h=1:size(eta0,2)
                Sleep_update=Sleep;
                for i=1:size(Waking_Time,2)
                    if delta_waking(i)>0
                        Sleep_update(Waking_Time(i)+1-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))):Waking_Time(i)+1,2)=0;
                    elseif delta_waking(i)<0
                        Sleep_update(Waking_Time(i):Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),2)=1;
                    end
                end
                for i=1:size(Sleeping_Time,2)
                    if delta_sleeping(i)<0
                        Sleep_update(Sleeping_Time(i):Sleeping_Time(i)-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                    elseif delta_sleeping(i)>0
                        Sleep_update(Sleeping_Time(i)+1-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))):Sleeping_Time(i)+1,2)=1;
                    end
                end
                sim('Gradient_Descent_Sleep_Advance_Line_Search.slx')
                T(h)=theta_plus(end,1);
                error(h)=Error(end,2);
            end
            [M,N]=min(T);
            Num=find(T==T(N));
            [MM,NN]=min(error(Num));
            eta=eta0(Num(NN));
            for i=1:size(Waking_Time,2)
                if delta_waking(i)>0
                    Sleep(Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))):Waking_Time(i),2)=0;
                elseif delta_waking(i)<0
                    Sleep(Waking_Time(i):Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))),2)=1;
                end
            end
            for i=1:size(Sleeping_Time,2)
                if delta_sleeping(i)<0
                    Sleep(Sleeping_Time(i):Sleeping_Time(i)-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                elseif delta_sleeping(i)>0
                    Sleep(Sleeping_Time(i)+1-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))):Sleeping_Time(i)+1,2)=1;
                end
            end
            sim('Gradient_Descent_Sleep_Advance.slx')
            T_opt=theta(end,1)
        else
            for h=1:size(eta0,2)
                Sleep_update=Sleep;
                for i=1:size(Waking_Time,2)
                    if delta_waking(i)>0
                        Sleep_update(max([Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),1]):Waking_Time(i),2)=0;
                    elseif delta_waking(i)<0
                        Sleep_update(Waking_Time(i):Waking_Time(i)-round(eta0(h)*delta_waking(i)/abs(delta_waking(i))),2)=1;
                    end
                end
                for i=1:size(Sleeping_Time,2)
                    if delta_sleeping(i)<0
                        Sleep_update(Sleeping_Time(i):Sleeping_Time(i)-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                    elseif delta_sleeping(i)>0
                        Sleep_update(max([Sleeping_Time(i)+1-round(eta0(h)*delta_sleeping(i)/abs(delta_sleeping(i))),1]):Sleeping_Time(i)+1,2)=1;
                    end
                end
                sim('Gradient_Descent_Sleep_Delay_Line_Search.slx')
                T(h)=theta_plus(end,1);
                error(h)=Error(end,2);
            end
            [M,N]=min(T);
            Num=find(T==T(N));
            [MM,NN]=min(error(Num));
            eta=eta0(Num(NN));
            for i=1:size(Waking_Time,2)
                if delta_waking(i)>0
                    Sleep(Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))):Waking_Time(i),2)=0;
                elseif delta_waking(i)<0
                    Sleep(Waking_Time(i):Waking_Time(i)-round(eta*delta_waking(i)/abs(delta_waking(i))),2)=1;
                end
            end
            for i=1:size(Sleeping_Time,2)
                if delta_sleeping(i)<0
                    Sleep(Sleeping_Time(i):Sleeping_Time(i)-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))),2)=0;
                elseif delta_sleeping(i)>0
                    Sleep(Sleeping_Time(i)+1-round(eta*delta_sleeping(i)/abs(delta_sleeping(i))):Sleeping_Time(i)+1,2)=1;
                end
            end
            sim('Gradient_Descent_Sleep_Delay.slx')
            T_opt=theta(end,1)
        end
        %  clear theta u H
        T_optimal(j+1)=T_opt;
        Sleep=sleep;
        if eta==0
            break
        end
    end
    N=min([size(R1,1),size(theta,1)]);
    filename=strcat(num2str(n),'_shift_Controllable_Sleep_9am_1000lux','.mat');
    save(filename,'sleep','theta','u','H','R1','R2','T_optimal')
    T_Controllable(n)=T_optimal(end);
    figure (n)
    subplot(2,1,1)
    plot(Periodic_Solution(:,1),Periodic_Solution(:,2)/2/pi,'b','linewidth',2)
    hold on
    plot(Periodic_Solution(:,1),Periodic_Solution(:,3),'b--','linewidth',2)
    hold on
    plot(theta(:,1),theta(:,2)/2/pi,'k','linewidth',2)
    hold on
    plot(H(:,1),H(:,2),'k--','linewidth',2)
    hold on
    plot(u(:,1),u(:,2),'r','linewidth',2)
    grid on
    xlabel('time(hours)');
    axis([0 theta(end,1) 0 1.2])
    legend('\theta_{ref}','H_{ref}','\theta','H','u')
    set(gca,'fontsize',20)
    subplot(2,1,2)
    plot(R1(1:N,1),R1(1:N,2).*(1-sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20,'b','linewidth',2)
    hold on
    plot(u(:,1),u(:,2).*(1-sleep(:,2)),'r','linewidth',2)
    grid on
    xlabel('time(hours)');
    axis([0 theta(end,1) min(R1(1:N,2).*(1-sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20) max(R1(1:N,2).*(1-Sleep(1:N,2)).*interp1q(f(:,1),f(:,2),theta(1:N,2))/20)])
    legend('P_1(1-\beta)f(\theta)','u')
    set(gca,'fontsize',20)
end
