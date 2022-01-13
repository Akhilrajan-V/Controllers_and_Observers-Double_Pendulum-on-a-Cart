%% LQG for Linear System
 

%%System Model
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

B = [0; 1/M; 0; 1/M*l1; 0; 1/M*l2];
C = [1 0 1 0 1 0];

C1 = [1 0 0 0 0 0;
     0 0 1 0 0 0;
     0 0 0 0 1 0;];
    

D = 0;
QXU = 0.1*eye(7);
QWV = eye(7);
 
sys = ss(A, B,C, D);


KLQG = lqg(sys,QXU,QWV);
 
sys1 = ss(KLQG.A, KLQG.B, C1, 0);
 
X0 =[0;
    0;
    3.141;
    0;
    0;
    0;];
 
 figure
 initial(sys1,X0)%LQG initial condition
 
 figure 
 step(sys1)%LQG step input response 
 
 