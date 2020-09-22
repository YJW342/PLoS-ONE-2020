% clear all
% close all
% clc
% tol=0.01;
% Theta=[0:0.5:2*pi];
% Phi=[0:0.2:pi];
% T=[];error=[];
% for i=1:50:2400
%     load('Periodic_Solution_JFK_I_1000lux.mat');
%     Periodic_Solution=Periodic_Solution(i:end,:);
%     Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);
%     for h=1:size(Theta,2)
%         for k=1:size(Phi,2)
%             x0=Periodic_Solution(1,2)+sqrt(tol)*sin(Phi(k))*cos(Theta(h));
%             xc0=Periodic_Solution(1,3)+sqrt(tol)*sin(Phi(k))*sin(Theta(h));
%             H0=Periodic_Solution(1,4)+sqrt(tol)*cos(Phi(k));
%             n0=rand;
%             Sleep0=Periodic_Solution(1,6);
%             sim('Open_Loop_shift_notol_initial.mdl');
%             [time,portion]=terminal(Error,tol);
%             T=[T;[time,portion,max(Error(:,2))]];
%             error=[error;Error2];
%             hold on
%             plot(Error(:,1),Error(:,2))
%             grid on
%         end
%     end
%     filename=strcat('Error_',num2str(i),'.fig');
%     savefig(filename)
%     filename=strcat('Error_',num2str(i),'.mat');
%     save(filename,'T','error');
%     close all
%     T=[];
%     error=[];
% end
% 
% TT=[];
% for i=1:50:2400
%    filename=strcat('Error_',num2str(i),'.mat');
%    load(filename);
%    TT=[TT;T];
% end
% histogram(TT);


% absolute uniformly
clear all
close all
clc
tol=0.01;
Theta=[0:0.5:2*pi];
Phi=[0:0.2:pi];
T=[];error=[];
for i=1:1000
    N=round(2400*rand);
    if N==0 
        N=1;
    end
    load('Periodic_Solution_JFK_I_1000lux.mat');
    Periodic_Solution=Periodic_Solution(N:end,:);
    Periodic_Solution(:,1)=Periodic_Solution(:,1)-Periodic_Solution(1,1);    
    Theta=2*pi*rand;
    Phi=pi*rand;
    x0=Periodic_Solution(1,2)+sqrt(tol)*sin(Phi)*cos(Theta);
    xc0=Periodic_Solution(1,3)+sqrt(tol)*sin(Phi)*sin(Theta);
    H0=Periodic_Solution(1,4)+sqrt(tol)*cos(Phi);
    n0=rand;
    Sleep0=Periodic_Solution(1,6);
    sim('Open_Loop_shift_notol_initial.mdl');
    %[time,portion]=terminal(Error,tol);
    T=[T;[max(Error(:,2))]];
    hold on
    plot(Error(:,1),Error(:,2))
    grid on
end
histogram(T);

