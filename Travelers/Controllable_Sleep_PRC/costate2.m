function dR2 = costate2(t,R2)
global Sleep u theta
Sleep_p=interp1q(Sleep(:,1),Sleep(:,2),t);
u_p=interp1q(u(:,1),u(:,2),t);
theta_p=interp1q(theta(:,1),theta(:,2),t);
dR2=R2*[(1-Sleep_p)/18.2+Sleep_p/4.2];