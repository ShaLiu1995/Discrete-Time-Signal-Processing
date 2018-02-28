% A
w1 = randn(1, 1024);
w2 = sqrt(1/32) * randn(1, 1024);
h = (1/8) * ones(1, 8);
NFFT = 128;

H = fftshift(fft(h, NFFT));
magHdb = 20*log10(magH);
magH = abs(H);
phaH = angle(H)*180/pi;

figure;
plot(f, magHdb); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
title('Magnitude of H(k)');
xlabel('frequency (Hz)');  
ylabel('Magnitude (dB)');

figure;
plot(f, magH); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
title('Magnitude of H(k)');
xlabel('frequency (Hz)');  
ylabel('Magnitude');

figure;
plot(f, phaH); grid;
ylim([-180 180]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
title('Phase of H(k)');
xlabel('frequency (Hz)');  
ylabel('Phase (degrees)');

% B
beta = 5;
window = kaiser(NFFT, beta);
% overlap = 64;
fs = 1;
M = 128;
U = sum(window.^2) / M;
f = -0.5*fs : fs/NFFT : (0.5-1/NFFT)*fs;

% [Sw1w1,f] = pwelch(w1, window, overlap, NFFT, fs, 'centered');
W1 = zeros(15, 128);
lo = 1;
hi = 128;
for i = 1:15
   W1(i,:) = fftshift(fft(w1(lo:hi).*window', NFFT));
   lo = lo + 64;
   hi = hi + 64;
end

Sw1w1 = sum(abs(W1).^2, 1)./15;
Sw1w1 = Sw1w1/(fs*M*U);

figure;
plot(f, 10*log10(Sw1w1)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Sw1w1(f)');

figure;
plot(f, Sw1w1); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Sw1w1(f)');

% [Sw2w2,f] = pwelch(w2, window, overlap, NFFT, fs, 'centered');

W2 = zeros(15, 128);
lo = 1;
hi = 128;
for i = 1:15
   W2(i,:) = fftshift(fft(w2(lo:hi).*window', NFFT));
   lo = lo + 64;
   hi = hi + 64;
end

Sw2w2 = sum(abs(W2).^2, 1)./15;
Sw2w2 = Sw2w2/(fs*M*U);

figure;
plot(f, 10*log10(Sw2w2)); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Sw2w2(f)');

figure;
plot(f, Sw2w2); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Sw2w2(f)');

% Syy
y = filter(h,1,w1);
% [Syy,f] = pwelch(y, window, overlap, NFFT, fs, 'centered');

Y = zeros(15, 128);
lo = 1;
hi = 128;
for i = 1:15
   Y(i,:) = fftshift(fft(y(lo:hi).*window', NFFT));
   lo = lo + 64;
   hi = hi + 64;
end

Syy = sum(abs(Y).^2, 1)./15;
Syy = Syy/(fs*M*U);

figure;
plot(f, 10*log10(Syy)); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Syy(f)');

figure;
plot(f, Syy); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Syy(f)');

r = y + w2;
% [Srr,f] = pwelch(r, window, overlap, NFFT, fs, 'centered');

R = zeros(15, 128);
lo = 1;
hi = 128;
for i = 1:15
   R(i,:) = fftshift(fft(r(lo:hi).*window', NFFT));
   lo = lo + 64;
   hi = hi + 64;
end

Srr = sum(abs(R).^2, 1)./15;
Srr = Srr/(fs*M*U);

figure;
plot(f, 10*log10(Srr)); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Srr(f)');

figure;
plot(f, Srr); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Srr(f)');

% [Syw1,f] = cpsd(y, w1, [], [], NFFT, fs, 'centered');
Syw1 = sum(Y.*conj(W1), 1)./15;
Syw1 = Syw1/(fs*M*U);

figure;
plot(f, 10*log10(abs(Syw1))); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Syw1(f)');

figure;
plot(f, abs(Syw1)); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Syw1(f)');

figure;
plot(f, angle(Syw1)*180/pi); grid;
ylim([-180 180]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Phase (degrees)');
title('Power spectral estimates Syw1(f)');


% [Srw1,f] = cpsd(r, w1, [], [], NFFT, fs, 'centered');
Srw1 = sum(R.*conj(W1), 1)./15;
Srw1 = Srw1/(fs*M*U);

figure;
plot(f, 10*log10(abs(Srw1))); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power spectral estimates Srw1(f)');

figure;
plot(f, abs(Srw1)); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Power spectral estimates Srw1(f)');

figure;
plot(f, angle(Srw1)*180/pi); grid;
ylim([-180 180]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Phase (degrees)');
title('Power spectral estimates Srw1(f)');

% D
x = [0 0.1875 0.3125];
Hw1y = Syw1./Sw1w1;

figure;
plot(f, 10*log10(abs(Hw1y))); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Transfer function estimates Hw1y(f)');
hold on;
y=[10*log10(abs(Hw1y(65))) 10*log10(abs(Hw1y(89))) 10*log10(abs(Hw1y(105)))];
yneg = [1.3 1.3 1.3];
ypos = [1.1 1.1 1.1];
errorbar(x,y,yneg,ypos,'.');
hold off;

figure;
plot(f, abs(Hw1y)); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Transfer function estimates Hw1y(f)');

figure;
plot(f, angle(Hw1y)*180/pi); grid;
ylim([-180 180]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Phase (degrees)');
title('Transfer function estimates Hw1y(f)');
hold on
y=[angle(Hw1y(65))./pi.*180 angle(Hw1y(89))./pi.*180 angle(Hw1y(105))./pi.*180];
yneg = [8 8 8];
ypos = [8 8 8];
errorbar(x,y,yneg,ypos,'.');
hold off;

Hw1r = Srw1./Sw1w1;

figure;
plot(f, 10*log10(abs(Hw1r))); grid;
ylim([-30 5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Transfer function estimates Hw1r(f)');
hold on;
y = [10*log10(abs(Hw1r(65))) 10*log10(abs(Hw1r(89))) 10*log10(abs(Hw1r(105)))];
yneg = [1.3 2 4.5];
ypos = [1.1 1.6 3];
errorbar(x,y,yneg,ypos,'.');
hold off;

figure;
plot(f, abs(Hw1r)); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Transfer function estimates Hw1r(f)');

figure;
plot(f, angle(Hw1r)*180/pi); grid;
ylim([-180 180]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Phase (degrees)');
title('Transfer function estimates Hw1r(f)');
hold on
y=[angle(Hw1r(65))./pi.*180 angle(Hw1r(89))./pi.*180 angle(Hw1r(105))./pi.*180];
yneg = [8 12 24];
ypos = [8 12 24];
errorbar(x,y,yneg,ypos,'.');
hold off;

% E
Gamma2w1y=((abs(Syw1)).^2)./(Syy.*Sw1w1);
figure;
plot(f, Gamma2w1y); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Magnitude-squared coherence function estimates Gamma^2w1y(f)');
hold on
y=[0.9 0.9 0.9];
yneg = [0.09 0.09 0.09];
ypos = [0.04 0.04 0.04];
errorbar(x,y,yneg,ypos,'.');
hold off;

Gamma2w1r=((abs(Srw1)).^2)./(Srr.*Sw1w1);
figure;
plot(f, Gamma2w1r); grid;
ylim([0 1.5]);
set(gca,'xtick', [-0.5 -0.3125 -0.1875 0 0.1875 0.3125 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude');
title('Magnitude-squared coherence function estimates Gamma^2w1r(f)');
hold on
y = [0.9 0.8 0.5];
yneg = [0.09 0.15 0.25];
ypos = [0.04 0.08 0.17];
errorbar(x,y,yneg,ypos,'.');
hold off;
