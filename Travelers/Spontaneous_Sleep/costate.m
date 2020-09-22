function dR = costate(t,R)
global Sleep u x xc I
Ac=0.1333;tr=18.2;td=4.2;omega0=2*pi/24.2;beta=0.0075;I0=9500;
q = 1/3;taux = 24.2;k = 0.55;mu=0.13;G=33.75;alpha0=0.05;p=0.5;
Sleep_p=interp1q(Sleep(:,1),Sleep(:,2),t);
u_p=interp1q(u(:,1),u(:,2),t);
x_p=interp1q(x(:,1),x(:,2),t);
xc_p=interp1q(xc(:,1),xc(:,2),t);
I_p=interp1q(I(:,1),I(:,2),t);
Alpha_p=alpha0*(I_p/I0)^p;
% H_p=interp1q(H(:,1),H(:,2),t);
dR=[-pi*R(1)/12*(mu*(1/3+4*x_p^2-256*x_p^6/15)-0.4*(1-0.4*xc_p)*(1-Sleep_p)*u_p)-pi*R(2)/12*(-0.4*q*xc_p*(1-Sleep_p)*(1-0.4*xc_p)*u_p-(24/0.99729/taux)^2-k*(1-0.8*x_p)*(1-0.4*xc_p)*(1-Sleep_p)*u_p);
    -pi*R(1)/12*(1-0.4*(1-0.4*x_p)*(1-Sleep_p)*u_p)-pi*R(2)/12*(q*(1-0.8*xc_p)*(1-0.4*x_p)*(1-Sleep_p)*u_p+0.4*k*x_p*(1-0.4*x_p)*(1-Sleep_p)*u_p);
    R(3)*((1-Sleep_p)/tr+Sleep_p/td);
    pi*G*(1-0.4*x_p)*(1-0.4*xc_p)*(1-Sleep_p)*Alpha_p*(R(1)+R(2)*(q*xc_p-k*x_p))/12+60*R(4)*(Alpha_p*(1-Sleep_p)+beta)];