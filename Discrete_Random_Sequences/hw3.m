% I.A I.B
uniform = rand(1024,1);
gaussian = randn(1024,1);
figure;
h1 = histogram(uniform);
xlabel('value');
ylabel('Number of data points');
title ('Histogram of uniform distribution');
xlim([-0.2 1.2]);
grid;
figure;
h2 = histogram(gaussian);
xlabel('value');
ylabel('Number of data points');
title ('Histogram of Gaussian distribution');
grid;

% I.C
uniform_256 = uniform(1:256);
gaussian_256 = gaussian(1:256);
[cxx_u, lag_u] = xcorr(uniform_256, 15, 'biased');
figure;
plot(lag_u, cxx_u);
grid;
xlabel('corelation lag m');
ylabel('cxx(m)');
title ('Autocorrelation of uniform sequence');
xlim([0 15]);

[cxx_g,lag_g] = xcorr(gaussian_256, 15, 'biased');
figure;
plot(lag_g, cxx_g);
grid;
xlabel('corelation lag m');
ylabel('cxx(m)');
title ('Autocorrelation of Gaussian sequence');
xlim([0 15]);


% II.A.1
h = ones(1,8);
figure;
zplane(h);
title ('Zero-pole plot of the FIR filter');
grid;
% II.A.2
NFFT = 256;
h = [h, zeros(1, NFFT - 8)];
f = -0.5 : 1/NFFT : 0.5-1/NFFT;
H = fft(h,N);
H = fftshift(H);
figure;
plot(f,20*log10(abs(H)));     
title('FFT of h(n), NFFT = 256');
xlabel('f  (cycles/sample)');  
ylabel('Amplitude of H(k)  (dB)');
grid;

% II.B.1
figure;
plot(gaussian_256);
xlim([0 256]);
title('Gaussian sequence input');
xlabel('n');
ylabel('value');
grid;

% II.B.2
out_gaussian = filter(h, 1, gaussian_256);
figure;
plot(out_gaussian);
xlim([0 256]);
title('Gaussian sequence outpurput');
xlabel('n');
ylabel('value');
grid;

[cyy,lag_yy] = xcorr(out_gaussian, 15, 'biased');
figure;
plot(lag_yy, cyy);
grid;
xlabel('corelation lag m');
ylabel('cyy(m)');
title ('Cyy(m) of the output sequence');
xlim([0 15]);

[cxy,lag_xy] = xcorr(out_gaussian, gaussian_256, 15, 'biased');
figure;
plot(lag_xy, cxy);
grid;
xlabel('corelation lag m');
ylabel('Cxy(m)');
title ('Cxy(m) of the input and output sequence');

