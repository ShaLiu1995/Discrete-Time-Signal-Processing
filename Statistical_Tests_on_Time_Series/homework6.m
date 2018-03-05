% A
N = 1024;
x1 = randn(1, 1024);
n = 0: 1023;
A2 = sqrt(2);
x2 = x1 + A2*sin(pi/16*n);
A3 = sqrt(2*10^0.9);
x3 = x1 + A3*sin(pi/16*n);
x4 = [randn(1, 512) sqrt(4)*randn(1, 512)];

% B.1
figure; 
plot(n,x1); grid;
title('x1(n), mean=0, variance=1');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 1023]);

figure; 
plot(n,x2); grid;
title('x2(n), mean=0, variance=2');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 1023]);

figure; 
plot(n,x3); grid;
title('x3(n) mean=0, variance=9');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 1023]);

figure; 
plot(n,x4); grid;
title('x4(n), mean=0, variance=1');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 1023]);

% B.2
figure;
h1 = histfit(x1);
xlabel('value');
ylabel('Number of points');
title ('Probability density function estimate of x1');
grid;

figure;
h2 = histfit(x2);
xlabel('value');
ylabel('Number of points');
title ('Probability density function estimate of x2');
grid;

figure;
h3 = histfit(x3);
xlabel('value');
ylabel('Number of points');
title ('Probability density function estimate of x3');
grid;

figure;
h4 = histfit(x4);
xlabel('value');
ylabel('Number of points');
title ('Probability density function estimate of x4');
grid;

% B.3
fs = 1;
NFFT = 128;
overlap = 64;
beta = 5;
window = kaiser(NFFT, beta);

[X1,f] = pwelch(x1, window, overlap, NFFT, fs, 'centered');
figure;
plot(f, 10*log10(X1)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates x1');

[X2,f] = pwelch(x2, window, overlap, NFFT, fs, 'centered');
figure;
plot(f, 10*log10(X2)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates x2');

[X3,f] = pwelch(x3, window, overlap, NFFT, fs, 'centered');
figure;
plot(f, 10*log10(X3)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates x3');

[X4,f] = pwelch(x4, window, overlap, NFFT, fs, 'centered');
figure;
plot(f, 10*log10(X4)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates x4');


% C.1
x1 = x1(1:1000);
mean_sq_est_1 = zeros(1, 100);
for i = 0:99
    sq_sum = 0;
    for m = 1:10
        sq_sum = sq_sum + x1(10*i + m)^2;
    end
    mean_sq_est_1(i+1) = sq_sum/10;
end
i = 0:99;
median1 = median(mean_sq_est_1);
figure;
scatter(i, mean_sq_est_1); grid;
xlabel('i');
ylabel('value');
title('Sample mean-square estimate of x1');
hold on;
plot(i, median1*ones(1,100),'r--');
hold off;

x2 = x2(1:1000);
mean_sq_est_2 = zeros(1, 100);
for i = 0:99
    sq_sum = 0;
    for m = 1:10
        sq_sum = sq_sum + x2(10*i + m)^2;
    end
    mean_sq_est_2(i+1) = sq_sum/10;
end
i = 0:99;
median2 = median(mean_sq_est_2);
figure;
scatter(i, mean_sq_est_2); grid;
xlabel('i');
ylabel('value');
title('Sample mean-square estimate of x2');
hold on;
plot(i, median2*ones(1,100),'r--');
hold off;

x3 = x3(1:1000);
mean_sq_est_3 = zeros(1, 100);
for i = 0:99
    sq_sum = 0;
    for m = 1:10
        sq_sum = sq_sum + x3(10*i + m)^2;
    end
    mean_sq_est_3(i+1) = sq_sum/10;
end
i = 0:99;
median3 = median(mean_sq_est_3);
figure;
scatter(i, mean_sq_est_3); grid;
xlabel('i');
ylabel('value');
title('Sample mean-square estimate of x2');
hold on;
plot(i, median3*ones(1,100),'r--');
hold off;

x4 = x4(1:1000);
mean_sq_est_4 = zeros(1, 100);
for i = 0:99
    sq_sum = 0;
    for m = 1:10
        sq_sum = sq_sum + x4(10*i + m)^2;
    end
    mean_sq_est_4(i+1) = sq_sum/10;
end
i = 0:99;
median4 = median(mean_sq_est_4);
figure;
scatter(i, mean_sq_est_4); grid;
xlabel('i');
ylabel('value');
title('Sample mean-square estimate of x4');
hold on;
plot(i, median4*ones(1,100),'r--');
hold off;

count1 = 1;
for i = 1:99
    if ((mean_sq_est_1(i) - median1) * (mean_sq_est_1(i+1) - median1) < 0)
        count1 = count1 + 1;
    end
end
    
count2 = 1;
for i = 1:99
    if ((mean_sq_est_2(i) - median2) * (mean_sq_est_2(i+1) - median2) < 0)
        count2 = count2 + 1;
    end
end

count3 = 1;
for i = 1:99
    if ((mean_sq_est_3(i) - median3) * (mean_sq_est_3(i+1) - median3) < 0)
        count3 = count3 + 1;
    end
end

count4 = 1;
for i = 1:99
    if ((mean_sq_est_4(i) - median4) * (mean_sq_est_4(i+1) - median4) < 0)
        count4 = count4 + 1;
    end
end

% C.2
N = 1000;
K = 30;
n = K-3;
F = N/K*ones(1, K);
alpha = (29:-1:0)/30;
z_alpha = icdf('normal',1-alpha,0, 1);

mean1 = mean(x1);
std1 = std(x1);
chi1 = mean1 + std1*z_alpha;
f1 = zeros(1, K);
for i = 2:K
    f1(i) = length(find(x1>chi1(i-1) & x1<chi1(i)));
end
f1(1) = N - sum(f1);
abs_ff_1 = abs(f1 - F);
chi_sq_1 = abs_ff_1.^2./F;
chi_sq_1_sum = sum(chi_sq_1)

mean2 = mean(x2);
std2 = std(x2);
chi2 = mean2 + std2*z_alpha;
f2 = zeros(1, K);
for i = 2:K
    f2(i) = length(find(x2>=chi2(i-1) & x2<chi2(i)));
end
f2(1) = N - sum(f2);
abs_ff_2 = abs(f2 - F);
chi_sq_2 = abs_ff_2.^2./F;
chi_sq_2_sum = sum(chi_sq_2)

index = 1:30;
index = index';
alpha = alpha';
z_alpha = z_alpha';
chi1 = chi1';
chi2 = chi2';
f1 = f1';
f2 = f2';
F = F';
abs_ff_1 = abs_ff_1';
abs_ff_2 = abs_ff_2';
chi_sq_1 = chi_sq_1';
chi_sq_2 = chi_sq_2';
T = table(index,alpha,z_alpha,chi1,f1,F,abs_ff_1,chi_sq_1)
T = table(index,alpha,z_alpha,chi2,f2,F,abs_ff_2,chi_sq_2)



