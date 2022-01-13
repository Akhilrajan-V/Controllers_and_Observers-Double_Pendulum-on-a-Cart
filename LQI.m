%% LQ Integral for Constant reference tracking
% Applied on Linearized system 

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

Q = [1000 0 0 0 0 0;
     0 1000 0 0 0 0;
     0 0 1000 0 0 0;
     0 0 0 1000 0 0;
     0 0 0 0 100 0;
     0 0 0 0 0 1000;];
 
 R = 0.01;
 
 
 C = [1 0 0 0 0 0;];
    
     
 D = 0; 

 %%Initial Condition
 X0 =[0;
     0;
     3.141;
     0;
     0;
     0;];

k= lqr(A,B,Q,R);
kx =k(1:3);
ki =k(4:6);
 X_ref =[2;
        0;
        0;
        0;
        0;
        0;];
     
N_bar = -inv(C*inv(A-B*k)*B);
BB = N_bar*B;
sys = ss(A-B*k,BB,C,D);
opt = stepDataOptions('InputOffset',0,'StepAmplitude',2);
tFinal = 300;
step(sys,opt)
title('X Tracking a Constant Reference');

