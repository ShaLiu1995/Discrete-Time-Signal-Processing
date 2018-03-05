% I.A
% Generate the 4096-point signal
n = 0:4095;
Fs = 1000;
A1 = 100; A2 = 10;  A3 = 1;
f1 = 160; f2 = 237; f3 = 240;
x = A1*cos(2*pi*(f1/Fs)*n) + A2*cos(2*pi*(f2/Fs)*n) + A3*cos(2*pi*(f3/Fs)*n);
figure; 
plot(n,x); grid;
title('x(n) (n=0,...,255)');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 255]);

% I.B
% FFT of the windowed signal
beta = 7.85;
NFFT = 256;
window = kaiser(NFFT, beta);
x_window = x(1:256).*window';
X = fftshift(fft(x_window, NFFT));
magX = abs(X);
magXdb = 20*log10(magX);
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
figure;
plot(f, magXdb); grid;
ylim([0 80]); xlim([-500 500]);
set(gca,'xtick', [-500 -240 -160 0 160 240 500]);
title('Magnitude of X(k)');
xlabel('frequency (Hz)');  
ylabel('Magnitude (dB)');

% apply the rectangle window
% window = rectwin(NFFT);
% x_window = x(1:256).*window';
% X = fftshift(fft(x_window, NFFT));
% magX = abs(X);
% magXdb = 20*log10(magX);
% f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
% figure;
% plot(f, magXdb); grid;
% % ylim([0 80]); xlim([-500 500]);
% set(gca,'xtick', [-500 -240 -160 0 160 240 500]);
% title('Magnitude of X(k) (rectangle window)');
% xlabel('frequency (Hz)');  
% ylabel('Magnitude (dB)');

% II.A
% Design Low-pass filter
h = firpm(63,[0 40 85 500]/500,[1 1 0 0],[1 150]);
% fvtool(h ,1);
figure; 
plot(h); grid;
title('Amplitude of h(n)');
xlabel('n');  
ylabel('Amplitude');
xlim([0 63]);

% II.B
% FFT of the low-pass filter
NFFT = 1024;
h_zeros = [h zeros(1, NFFT-64)];
h_window = h_zeros.*rectwin(NFFT)';
H = fftshift(fft(h_window, NFFT));
magH = abs(H);
magHdb = 20*log10(magH);
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
figure;
plot(f, magHdb); grid;
set(gca,'xtick', [-500 -250 0 250 500]);
title('Magnitude of H(k)');
xlabel('Frequency (Hz)');  
ylabel('Magnitude (dB)');

% Blow up in 0<=|f|<=40Hz
figure;
plot(f, magHdb); grid;
xlim([-40 40]);
title('Magnitude of H(k) (0<=|f|<=40Hz)');
xlabel('Frequency (Hz)');  
ylabel('Magnitude (dB)');

% III.A
NFFT = 256;
f0 = 250;
x_exp = x.*exp(-1i*2*pi*(f0/Fs)*n);
x_exp_window = x_exp(1:256).*kaiser(NFFT, beta)';
X = fftshift(fft(x_exp_window, NFFT));
magX = abs(X);
magXdb = 20*log10(magX);
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
figure;
plot(f, magXdb); grid;
ylim([0 80]); xlim([-500 500]);
% set(gca,'xtick', [-500 -250 0 250 500]);
set(gca,'xtick', [ -490 -410 -250 -90 -10 250 500]);
title('Magnitude of X(k) after multiplying exp(-jw_0n)');
xlabel('Frequency (Hz)');  
ylabel('Magnitude (dB)');

% III.B
NFFT = 256;
x_lpf = filter(h, 1, x_exp);
x_lpf_window = x_lpf(256:511).*kaiser(NFFT, beta)';
X = fftshift(fft(x_lpf_window, NFFT));
magX = abs(X);
magXdb = 20*log10(magX);
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
figure;
plot(f, magXdb); grid;
ylim([0 80]); xlim([-500 500]);
% set(gca,'xtick', [-500 -250 0 250 500]);
set(gca,'xtick', [ -490 -410 -250 -90 -10 250 500]);
title('Magnitude of X(k) after filtering');
xlabel('Frequency (Hz)');  
ylabel('Magnitude (dB)');

% III.C
x_desample = downsample(x_lpf, 8);
figure;
plot(abs(x_desample)); grid;
title("x'(n) after desampling");
xlabel('n'); 
ylabel('Magnitude |x(n)|');
xlim([0 511]);

% IV
NFFT = 256;
FsPrime = Fs/8;
x_desample_window = x_desample(256:511).*kaiser(NFFT, beta)';
X = fftshift(fft(x_desample_window, NFFT));
magX = abs(X);
magXdb = 20*log10(magX);
f = -0.5*FsPrime : FsPrime/NFFT : (0.5-1/NFFT)*FsPrime;
figure;
plot(f, magXdb); grid;
ylim([-40 60]); 
xlim([-63 63]);
set(gca,'xtick', [-63 0 63]);
set(gca,'xtick', [-62 -35 -10 0 10 35 62]);
title("Magnitude of X'(k)  (LPF weight ratio 1:150)");
xlabel('Frequency (Hz)');  
ylabel('Magnitude (dB)');