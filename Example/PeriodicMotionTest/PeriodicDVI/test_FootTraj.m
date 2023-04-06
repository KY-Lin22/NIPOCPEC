
k_N = 0.5;
k_T = 0.25;

rho = 1: -0.01: -1;

pN = k_N * 1/pi * (cos(pi*rho) + 1);
pT = k_T * 2/pi * (sin(-pi/2*rho) + 1);

figure(1)
subplot(2,1,1)
plot(pN, 'k')
legend('posN')
xlabel('time')
subplot(2,1,2)
plot(pT, 'k')
legend('posT')
xlabel('time')
figure(2)
plot(pT, pN)
%%
clear all
clc
a1 = 0.2;
a2 = 50;
k_T = 0.25;

rho = 1: -0.01: -1;
num = size(rho, 2);
pN = zeros(1, num);
for n = 1 : num
    pN(n) = a1 / (1 + a2 * (rho(n))^2);
end

pT = k_T * 2/pi * (sin(-pi/2*rho) + 1);
plot(pT, pN)

