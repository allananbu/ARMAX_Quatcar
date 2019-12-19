clc
clear all
close all

%define system parameters
m1=466.5;
m2=49.8;
k1=5700;
k2=135000;
b1=290;
b2=1400;
k1=5.51;
k2=10;

%define the matrices for MIMO system
A=[-b1/m1 b1/m1 -k1/m1 k1/m1; 
    b1/m2 -(b1+b2)/m2 k1/m2 -(k1+k2)/m2; 
    1 0 0 0; 
    0 1 0 0];
B=[0 0;
    b2/m2 k2/m2;
    0 0;
    0 0];
C=[1 0 0 0];
D=zeros(1,2);
sys=ss(A,B,C,D);
t=0:0.1:10;
x0=[0.5;2.4;5;8];
u=zeros(2,101);
[y,x]=lsim(sys,u,t,x0);
plot(t,y);

%convert to transfer function
[num,den]=ss2tf(A,B,C,D,1);

sys=tf(num,den);
sys_d=c2d(sys,1,'zoh');

%build idpoly model 
%take A,B,C polynomials from the discrete model
num= [0 0.3939 -0.7711 0.36817 0.009080]; %A polynomial
den= [1 -2.5825 2.1755 -0.5929 9.87523215909130e-16]; %B polynomial
mod = idpoly(den,num); 
step(mod);
u=idinput(1000,'prbs',[0 0.5],[0 1]); %PRBS input
N=size(u);
y=sim(mod,u);

%build estimator model
order=[4 4 1 1];
data=[y u];
SYS=armax(data,order);
step(SYS);