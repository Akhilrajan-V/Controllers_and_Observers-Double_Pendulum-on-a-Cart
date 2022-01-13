%% Luenberger Observer for the Non Linear Model 

X0 =[0;
     0;
     3.141;
     0;
     0;
     0;];
time = 0:0.1:350;
[t1,y2] = ode45(@luenberger_nl,time,X0); % using ode45 function for Non Linear ODE

figure
plot(t1,y2(:,1))
title('Observation of X (Non Linear)')
xlabel('time(Seconds)')
ylabel('Distance(m)')
grid on

figure
plot(t1,y2(:,3))
title('Observation of theta1 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

figure
plot(t1,y2(:,5))
title('Observation of theta2 (Non Linear)')
xlabel('time(Seconds)')
ylabel('radians')
grid on

function dydt = luenberger_nl(t,X0)
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
C1 = [1 0 0 0 0 0; 0 0 1 0 0 0; 0 0 0 0 1 0]; % C matrix for X, theta1 and theat2
 
 
  Q=[1000 0 0 0 0 0;
     0 1000 0 0 0 0;
     0 0 1000 0 0 0;
     0 0 0 1000 0 0;
     0 0 0 0 100 0;
     0 0 0 0 0 1000];
 
  R = 0.01;
   
  k = lqr(A,B,Q,R);
  
  F=-k*X0;
        
  pole = place(A',C1',[-20 -21 -25 -24 -23 -26;]);

  L = pole';

  y1 = [X0(1); X0(3); X0(5)];

  Error_comp = L*(y1-C1*X0);
 
 dydt=zeros(6,1);
 % y(1)=x; y(2)=xdot; y(3)=theta1;   y(4)=theta1dot;  y(5)=theta2;    y(6)=theta2dot;
 dydt(1)= X0(2)+Error_comp(1); %XD
 dydt(2)= (F-m1*g*sin(X0(3))*cos(X0(3))-m1*l1*(X0(4)^2)*sin(X0(3))-m2*g*sin(X0(5))*cos(X0(5))-m2*l2*(X0(6)^2)*sin(X0(5)))/(M+m1*((sin(X0(3)))^2)+m2*((sin(X0(5)))^2))+Error_comp(2);
 dydt(3)= X0(4)+Error_comp(3); %theta 1D
 dydt(4)= (dydt(2)*cos(X0(3))-g*(sin(X0(3))))/l1+Error_comp(4); %theta 1 Ddot;
 dydt(5)= X0(6)+Error_comp(5); %theta 2D
 dydt(6)= (dydt(2)*cos(X0(5))-g*(sin(X0(5))))/l2+Error_comp(6); %theta 2Ddot;
end