% Part A & B
% Plot h(n)
n = 0:63;
h = firpm(63,[.1 .9],[1 1],'hilbert');
% figure;
% stem(n,h); 
% title('Impulse Response h(n) of the Hilbert Transformer');
% xlabel('n');  
% ylabel('h(n)');
% xlim([0 63]);
% grid;

% Plot H(k)
N = 1024;
f = -0.5 : 1/N : 0.5-1/N;
H = fft(h,N);
H = fftshift(H);
% figure;
% plot(f,20*log10(abs(H)));     
% title('Amplitude of H(k)');
% xlabel('f  (cycles/sample)');  
% ylabel('Amplitude of H(k)  (dB)');
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% grid;

% figure;
% plot(f,angle(H)); 
% axis([-0.5 0.5 -4 4]);
% title('Phase of H(k)');
% xlabel('f  (cycles/sample)');  
% ylabel('Phase of H(k) (rad)'); 
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% grid;


% blow up the phase of H(k) 
% figure,  
% plot(f,angle(H));        
% title('Phase of H(k)');
% xlabel(' f  (cycles/sample)');  
% ylabel('Phase of H(k) (rad)'); 
% axis([-0.01 0.01 -4 4]);
% grid;

% Part D
% low-pass filter
hlp = firpm(63,[0 .5 .6 1.0],[1 1 0 0]);
% figure;
% stem(n,hlp);
% title('Impulse response of the low-pass filter');
% xlabel('n');  
% ylabel('hlp(n)');
% xlim([0 63]);
% grid;

% Plot H_lp(k)
Hlp = fft(hlp,N);
Hlp = fftshift(Hlp);
% figure; 
% plot(f, 20*log10(abs(Hlp)));        
% title('Magnitude of Hlp(k)');
% xlabel('f  (cycles/sample)');  
% ylabel('Magnitude of Hlp(k) (dB)');
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% grid;

% figure,  
% plot(f,angle(Hlp));        
% title('Phase of Hlp(k)');
% xlabel('f  (cycles/sample)');  
% ylabel('Phase of Hlp(k) (rad)');
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% grid;

% generate xr(n)
n = 0:1023;
A1 = 10;
A2 = 10^0.5; 
A3 = 1;
f1 = 0.05;
f2 = 0.075; 
f3 = 0.10;
fc = 0.25;
xr = A1*sin(2*pi*f1*n) + A2*sin(2*pi*f2*n) + A3*sin(2*pi*f3*n);

NFFT = 256;
window = kaiser(NFFT,10);
% window = hamming(NFFT);
f_256 = -0.5 : 1/NFFT : 0.5-1/NFFT;

% Xr(k)
Xr = fftshift(fft(xr(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(Xr)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title('Magnitude of Xr(k)');
% xlabel('f (cycles/sample)');
% ylabel('Magnitude of Xr(k) (dB)');
% grid;

% Xr'(k)
xrPrime = conv(xr, hlp);
XrPrime = fftshift(fft(xrPrime(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(XrPrime)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("Magnitude of Xr'(k)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude of Xr'(k) (dB)");
% grid;

% Xi(k)
xi = conv(xr, h);
Xi = fftshift(fft(xi(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(Xi)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("Magnitude of Xi(k)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude of Xi(k) (dB)");
% grid;

xrcos = xrPrime(64:end).*cos(2*pi*fc*n);
Xrcos = fftshift(fft(xrcos(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(Xrcos)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("FFT of xr'(n)*cos(2*pi*fc*n)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude (dB)");
% grid;

xisin = xi(64:end).*sin(2*pi*fc*n);
Xisin = fftshift(fft(xisin(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(Xisin)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("FFT of xi(n)*sin(2*pi*fc*n)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude (dB)");
% grid;

srUSB = xrcos - xisin;
SrUSB = fftshift(fft(srUSB(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(SrUSB)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("Magnitude of Sr(k) (USB)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude of Sr(k) (USB) (dB)");
% grid;

srLSB = xrcos + xisin;
SrLSB = fftshift(fft(srLSB(257:512).*window', NFFT));
% figure;
% plot(f_256, 20*log10(abs(SrLSB)));
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title("Magnitude of Sr(k) (LSB)");
% xlabel('f (cycles/sample)');
% ylabel("Magnitude of Sr(k) (LSB) (dB)");
% grid;

xn = xrPrime + 1i*xi;
Xn = fftshift(fft(xn(257:512).*window', NFFT));
figure;
plot(f_256, 20*log10(abs(Xn)));
ylim([0 80]);
xlim([-0.5 0.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
title("Magnitude of X(k)");
xlabel('f (cycles/sample)');
ylabel("Magnitude of X(k) (dB)");
grid;

sn = xn(64:end) .* exp(1i*2*pi*fc*n);
Sn = fftshift(fft(sn(257:512).*window', NFFT));
figure;
plot(f_256, 20*log10(abs(Sn)));
ylim([0 80]);
xlim([-0.5 0.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
title("Magnitude of S(k)");
xlabel('f (cycles/sample)');
ylabel("Magnitude of S(k) (dB)");
grid;