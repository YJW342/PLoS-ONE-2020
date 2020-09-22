clear all;
close all;
clc;
bdclose('all');
load('Periodic_Solution_JFK_I_1000lux.mat');
Imax=10000;umax=0.2208;
tol=0.01;Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;alpha0=0.05;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;p=0.5;I0=9500;beta=0.0075;
Initial_Time=1;
Periodic_Solution=Periodic_Solution(Initial_Time*100+1:end,:);
Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
time_shift=[1:23];Sleep0=0;
for nn=1:size(time_shift,2)
    T_optimal=[];
    n0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,7),time_shift(nn));
    x0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,2),time_shift(nn));
    xc0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,3),time_shift(nn));
    H0=interp1q(Periodic_Solution(:,1),Periodic_Solution(:,4),time_shift(nn));
    cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Controllable_Sleep_PRC'
    filename=strcat(num2str(nn),'_shift_Controllable_Sleep_7am_10000lux','.mat'); 
    load(filename);
    I_update=u;I_update(:,2)=I_update(:,2)*Imax/umax;
    Sleep_update=sleep;
    cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Controllable_Sleep'
    sim('Integration_Forward_JFK_Simplification.mdl');
    T_1st(nn)=x(end,1);
    
    cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Controllable_Sleep_2nd';
    load(filename);
    I_update=u;I_update(:,2)=I_update(:,2)*Imax/umax;
    Sleep_update=sleep;
    cd 'C:\Users\Jiawei Yin\Dropbox\2020.02.26-2020.03.03\Travelers\Controllable_Sleep'
    sim('Integration_Forward_JFK_Simplification.mdl');
    T_2nd(nn)=x(end,1);   
    
    filename=strcat(num2str(nn),'_shift_Controllable_7am_10000lux','.mat'); 
    load(filename);
    T_OPT(nn)=T_optimal(end);
end
plot(T_1st)
hold on
plot(T_2nd)
hold on
plot(T_OPT);
grid on