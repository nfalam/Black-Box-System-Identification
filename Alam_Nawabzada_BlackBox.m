clc
clear all 
close all

N = 2000;
K = 50;
A = 2;              % Change A = 0.2 for Low SNR for Noise Profile, H(q)
u = A*dprbs(N, K);
[y,Ts] = project_blackbox_data(u,'52');

%% Data samples into Data Identification (70%) and Data Validation (30%)
Ts = 0.05;
datan = iddata(y,u',Ts);
data = dtrend(datan);
figure(2)
plot(data)
grid on
datai=data(1:1400);                                                          %Identification data, Nid = 700 
datav=data(1401:N);                                                          %Validation data, Nval = 300

%%
save data.mat 'u' 'y' 'Ts' 'datan' 'data' 'datai' 'datav'
%load data.mat
M = 400;
Nid = 1400;
Nval = 700;
g = etfe(datai,M,Nid,Ts);
figure(3)
bode(g);

%% Calculating and plotting Impulse Response and Impulse Weights to find out delay, nk
L = 50;
figure(4)
w = cra(datai,L);
us = ones(size(w));                                                         %generate unit-step input signal
figure(5)
ys = decresp(w,us);                                                         %recovered step response weights
xlim ([0 50]);
t = linspace (0, 50*Ts, 51);
title ('Step Response by Convolution of input with Impulse Weights');

nk = 13;                                                                    %nk from impulse arbitrary response
figure(6)
Sdet= hanktest(w,nk,0); 

%% OE Model
nb = 2;                                                                     %Numerator Order (Zeroes)
nf = 2;                                                                     %Denominator Order (Poles)
nk = 13;                                                                    %Delays
M_OE = oe(datai,[nb nf nk]);

figure(7)
resid(data, M_OE)
present(M_OE)
figure(8)
compare (data,M_OE,1)
%disp('CHECK HERE')

%% Chi-Squared Test OE Model
[res,R] = resid(data, M_OE);
res_OE = get(res,'y');
m = 25;
[Sr_OE,xr_OE] = chisq(m,res_OE,0)
[Sru_OE,xru_OE] = chisq(m,res_OE,u',0)

%% Noise Profile H(q)
figure(9);
autocorr(res_OE, 30)
figure(10);
parcorr(res_OE, 30)

n = data.y;                      
Hq_model = armax(n,[1,2])
%Hq_model = ar(y,1)
resid(n,Hq_model);
[res,R] = resid(n, Hq_model);
[S,x] = chisq(25,res)

%% Model input Parameters for BJ and PEM models
na = 1;
nb = 2;
nc = 1;
nd = 2;
nf = 2;
nk = 13;

%% Box Jenkins
M_BJ = bj(datai,[nb nc nd nf nk])
figure(11)
resid(data, M_BJ)
present(M_BJ)
figure(12)
compare (data,M_BJ,1)

%% Chi squared test for BJ
[res,R] = resid(data, M_BJ);
res_BJ = get(res,'y');
m = 25;
[Sr_BJ,xr_BJ] = chisq(m,res_BJ,0)
[Sru_BJ,xru_BJ] = chisq(m,res_BJ,u',0)


%% PEM
M_PEM = pem(datai,[na nb nc nd nf nk])
figure(13)
resid(data, M_PEM)
present(M_PEM)
figure(14)
compare (data,M_PEM,1)

%% Chi squared test for PEM
[res,R] = resid(data, M_PEM);
res_PEM = get(res,'y');
m = 25;
[Sr_PEM,xr_PEM] = chisq(m,res_PEM,0)
[Sru_PEM,xru_PEM] = chisq(m,res_PEM,u',0)

%% Step Response for OE, BJ and PEM 
G_OE = d2c (M_OE, 'zoh');
G_BJ = d2c (M_BJ, 'zoh');
G_PEM = d2c (M_PEM, 'zoh');
figure(15)
plot(step(G_OE)); 
hold on
plot(step(G_BJ));
hold on
plot(step(G_PEM));
hold off
legend('OE', 'BJ', 'PEM');
title('Step Response');
xlabel('Time (seconds)');
ylabel ('Amplitude');

%% Pole Zero Map
[A,B,C,D,F] = polydata(M_BJ);
Hq_cl = tf(C,D,Ts);
Gq_cl = tf(B,F,Ts);
numC = 2.2*[1 -0.75];
denC = [1 -0.25];
Gc = tf(numC,denC,Ts);

Gq = feedback(Gq_cl,Gc,+1);
Hq = (1+Gq*Gc)*Hq_cl;

figure(16)
pzmap(Hq);
figure(17)
pzmap(Gq); 

G_OE
G_BJ
G_PEM

