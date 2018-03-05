clear all;
clc;
%ECE251A HW2
%Generate a hilbert transformer
N = 1024;
NFFT = 256;

n = 0:63;
h = firpm(63,[0.1, 0.9],[1,1],'hilbert');
figure,
stem(n,h);
axis([0,63,-0.9,0.9]);
title('Figure 1: Impulse Response of the Hilbert Transformer');
xlabel('n');  
ylabel('h(n)');
Hk = fftshift(fft(h,N));
interval = 1/N;
f = -0.5:interval:0.5-interval;

% Plot magnitude and phase of Hk.  
figure,  
plot(f,20*log10(abs(Hk)));        
title('Figure 2: Magnitude of H(k)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of H(k)(dB)');
grid on;
figure,  
plot(f,angle(Hk));        
title('Figure 3: Phase of H(k)');
xlabel(' f(cycles/sample)');  
ylabel('Phase of H(k)'); grid on;

% blow up the phase plot in the vicinity of f = 0 cycles/sample  
figure,  
plot(f,20*log10(abs(Hk)));        
title('Figure 4: Magnitude of H(k)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of H(k)(dB)');
axis([-0.1,0.1,-30,5]);
grid on;
figure,  
plot(f,angle(Hk));        
title('Figure 5: Phase of H(k)');
xlabel(' f(cycles/sample)');  
ylabel('Phase of H(k)'); 
axis([-0.1,0.1,-4,4]);

%D
%generate a low pass filter
h1 = firpm(63,[0,0.8,0.9,1.0],[1,1,0,0]);
figure,
stem(n,h1);
axis([0,63,-0.9,0.9]);
title('Figure 6: Impulse Response of the Low Pass Filter');
xlabel('n');  
ylabel('h1(n)');
H1k = fftshift(fft(h1,N));

% Plot magnitude and phase of H1k.  
figure,  
plot(f,20*log10(abs(H1k)));        
title('Figure 7: Magnitude of H1(k)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of H1(k)(dB)');
grid on;
figure,  
plot(f,angle(H1k));        
title('Figure 8: Phase of H1(k)');
xlabel(' f(cycles/sample)');  
ylabel('Phase of H1(k)'); grid on;

%window
win = (blackman(1024))';
A1 = 10^(20/20);
A2 = 10^(10/20); 
A3 = 10^(0/20);
f1 = 0.05;
f2 = 0.075; 
f3 = 0.10;
fc = 0.25;

%xr(n)
n=0:1023;
xr=A1*sin(2*pi*f1*n)+A2*sin(2*pi*f2*n)+A3*sin(2*pi*f3*n);
%xr = xr(64:end);
Xrk = fftshift(fft(xr.*win,NFFT));
interval = 1/NFFT;
f = -0.5:interval:0.5-interval;
% Plot magnitude phase of Xrk.  
figure,  
plot(f,20*log10(abs(Xrk)));        
title('Figure 9: Magnitude of Xrk');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Hr(k)(dB)');
grid on;

%xi(n)
xi = conv(xr,h);
xi = xi(64:end);
Xik = fftshift(fft(xi.*win,NFFT));
% Plot magnitude phase of Xik.  
figure,  
plot(f,20*log10(abs(Xik)));        
title('Figure 10: Magnitude of Xik');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Xi(k)(dB)');
grid on;

%xr'(n)
xr1 = conv(xr,h1);
xr1 = xr1(64:end);
Xr1k = fftshift(fft(xr1.*win,NFFT));
% Plot magnitude of X'rk.  
figure,  
plot(f,20*log10(abs(Xr1k)));        
title('Figure 11: Magnitude of Xr1k');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Xr1(k)(dB)');
grid on;

%xr’(n)*cos(2?fcn)
xr2 = xr1.*cos(2*pi*fc*n);
Xr2k = fftshift(fft(xr2.*win,NFFT));
% Plot magnitude of Xr2k.  
figure,  
plot(f,20*log10(abs(Xr2k)));        
title('Figure 12: Magnitude of FFT of xr’(n)*cos(2?fcn)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Xr2k (dB)');
grid on;

%xi(n)*sin(2?fcn)
xi2 = xi.*sin(2*pi*fc*n);
Xi2k = fftshift(fft(xi2.*win,NFFT));
% Plot magnitude of Xi2k.  
figure,  
plot(f,20*log10(abs(Xi2k)));        
title('Figure 13: Magnitude of FFT of xi(n)*sin(2?fcn)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Xi2k (dB)');
grid on;

%sr(n) USB
sru = xr2-xi2;
Sruk = fftshift(fft(sru.*win,NFFT));
% Plot magnitude of Sruk.  
figure,  
plot(f,20*log10(abs(Sruk)));        
title('Figure 14: Magnitude of FFT of sr(n)(USB)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Sruk (dB)');
grid on;

%sr(n) LSB
srl = xr2+xi2;
Srlk = fftshift(fft(srl.*win,NFFT));
% Plot magnitude of Srlk.  
figure,  
plot(f,20*log10(abs(Srlk)));        
title('Figure 15: Magnitude of FFT of sr(n)(LSB)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of Srlk (dB)');
grid on;

%x(n)=x'r(n)+jxi(n)
x = xr1 + sqrt(-1)*xi;
X = fftshift(fft(x.*win,NFFT));
% Plot magnitude of X.  
figure,  
plot(f,20*log10(abs(X)));        
title('Figure 16: Magnitude of FFT of x(n)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of X (dB)');
grid on;

%s(n)=x(n)exp(+j2pifcn)
s = x.*exp(sqrt(-1)*2*pi*fc*n);
S = fftshift(fft(s.*win,NFFT));
% Plot magnitude of S.  
figure,  
plot(f,20*log10(abs(S)));        
title('Figure 17: Magnitude of FFT of s(n)');
xlabel(' f(cycles/sample)');  
ylabel('Magnitude of S (dB)');
grid on;