N = 1024;
x1 = randn(1, 1024);
n = 0: 1023;
A2 = sqrt(2);
x2 = x1 + A2*sin(pi/16*n);
A3 = sqrt(2*10^0.9);
x3 = x1 + A3*sin(pi/16*n);
x4 = [randn(1, 512) sqrt(4)*randn(1, 512)];

figure; 
plot(n,x1); grid;
title('x1(n), mean=0, variance=1');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 1023]);