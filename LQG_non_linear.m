%% LQG Controller for Non Linear Model

% Initial Condition
X0 =[0;
     0;
     3.141;
     0;
     0;
     0;];
time = 0:0.1:350;
[t1,y2] = ode45(@lqg_nl,time,X0); 

figure
plot(t1,y2(:,1))
title('Response of X (Non Linear)')
xlabel('time(Seconds)')
ylabel('Distance(m)')
grid on

figure
plot(t1,y2(:,3))
title('Response of theta1 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

figure
plot(t1,y2(:,5))
title('Response of theta2 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

function dydt = lqg_nl(t,y0)
M = 1000;
m1 = 100;
m2 = 100;
l1 = 20;
l2 = 10;
g = 9.8;

A = [0 1 0 0 0 0; 
    0 0 -(m1*g)/M 0 -(m2*g)/M 0;
    0 0 0 1 0 0;
    0 0 -(g/l1)-(g*m1)/M*l1 0 -(m2*g)/M*l1 0;
    0 0 0 0 0 1;
    0 0 -(m1*g)/M*l2 0 -(g/l2)-(g*m2)/M*l2 0];


B=[0; 1/M; 0; 1/(M*l1); 0; 1/(M*l2)];
C1 = [1 0 0 0 0 0]; % C matrix for X, theta1 and theat2

D = 0;    

 %%Initial Condition

 Q=[1000 0 0 0 0 0;
    0 1000 0 0 0 0;
    0 0 1000 0 0 0;
    0 0 0 1000 0 0;
    0 0 0 0 100 0;
    0 0 0 0 0 1000];
 
 R = 0.01;
 
  k = lqr(A,B,Q,R);
  F=-k*y0(1:6);
        
  pole = place(A',C1',[-20 -21 -25 -24 -23 -26;]);
  C = [1 0 0 0 0 0];

  L = pole';

  y1 = y0(1);

  QXU = 0.7*eye(6);
  QWV = 0.4;
 
  KLQG = lqr(A',C1',QXU,QWV);
   
 error_comp = KLQG'*(y1-C1*y0);
 
 dydt=zeros(6,1);
 % y(1)=x; y(2)=xdot; y(3)=theta1;   y(4)=theta1dot;  y(5)=theta2;    y(6)=theta2dot;
 dydt(1)= y0(2)+error_comp(1); %XD
 dydt(2)= (F-m1*g*sin(y0(3))*cos(y0(3))-m1*l1*(y0(4)^2)*sin(y0(3))-m2*g*sin(y0(5))*cos(y0(5))-m2*l2*(y0(6)^2)*sin(y0(5)))/(M+m1*((sin(y0(3)))^2)+m2*((sin(y0(5)))^2))+error_comp(2);
 dydt(3)= y0(4)+error_comp(3); %theta 1D
 dydt(4)= (dydt(2)*cos(y0(3))-g*(sin(y0(3))))/l1+error_comp(4); %theta 1 Ddot;
 dydt(5)= y0(6)+error_comp(5); %theta 2D
 dydt(6)= (dydt(2)*cos(y0(5))-g*(sin(y0(5))))/l2+error_comp(6); %theta 2Ddot;
end