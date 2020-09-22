function dR1 = costate1(t,R1)
global Sleep u f theta
Sleep_p=interp1q(Sleep(:,1),Sleep(:,2),t);
u_p=interp1q(u(:,1),u(:,2),t);
theta_p=interp1q(theta(:,1),theta(:,2),t);
dR1=-R1*interp1q(f(:,1),f(:,3),theta_p)*u_p*(1-Sleep_p);