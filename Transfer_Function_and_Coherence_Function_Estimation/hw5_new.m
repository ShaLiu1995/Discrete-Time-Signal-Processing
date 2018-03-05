% A
w1 = randn(1, 1024);
w2 = sqrt(1/32) * randn(1, 1024);
h = (1/8) * ones(1, 8);
NFFT = 128;
fs = 1;
f = -0.5 : 1/NFFT : (0.5-1/NFFT);
H = fftshift(fft(h, NFFT));
magH = abs(H);
magHdb = 20*log10(magH);
phaH = angle(H)*180/pi;

% figure;
% plot(f, magH); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title('Magnitude of H(k)');
% xlabel('frequency (Hz)');  
% ylabel('Magnitude');
% 
% figure;
% plot(f, magHdb); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title('Magnitude of H(k)');
% xlabel('frequency (Hz)');  
% ylabel('Magnitude (dB)');
% 
% figure;
% plot(f, phaH); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% title('Phase of H(k)');
% xlabel('frequency (Hz)');  
% ylabel('Phase (degrees)');

% B
beta = 5;
window = kaiser(NFFT, beta);
window = hamming(NFFT);
overlap = 64;
% fs = 1;
% M = 128;
% U = sum(window.^2) / M;
% segment the time series
% w1_reshape = reshape(w1,[128, 8])';
% % calculate the coherent average
% w1_avg = sum(w1_reshape, 1)./8;
% W1 = fft(w1_avg.*window' ,NFFT);
% W1 = fftshift(W1);
[Sw1w1,f] = pwelch(w1, window, overlap, NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(Sw1w1)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Sw1w1(f)');
% 
% figure;
% plot(f, Sw1w1); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Sw1w1(f)');

[Sw2w2,f] = pwelch(w2, window, overlap, NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(Sw2w2)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Sw2w2(f)');
% 
% figure;
% plot(f, Sw2w2); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Sw2w2(f)');

y = filter(h,1,w1);
[Syy,f] = pwelch(y, window, overlap, NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(Syy)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Syy(f)');
% 
% figure;
% plot(f, Syy); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Syy(f)');

r = y + w2;
[Srr,f] = pwelch(r, window, overlap, NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(Srr)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Srr(f)');
% 
% figure;
% plot(f, Srr); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Srr(f)');

[Syw1,f] = cpsd(y, w1, [], [], NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(abs(Syw1))); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Syw1(f)');
% 
% figure;
% plot(f, abs(Syw1)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Syw1(f)');

figure;
plot(f, angle(Syw1)*180/pi); grid;
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Phase (degrees)');
title('Power spectral estimates Syw1(f)');

[Srw1,f] = cpsd(r, w1, [], [], NFFT, fs, 'centered');

% figure;
% plot(f, 10*log10(abs(Srw1))); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude (dB)');
% title('Power spectral estimates Srw1(f)');
% 
% figure;
% plot(f, abs(Srw1)); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Magnitude');
% title('Power spectral estimates Srw1(f)');
% 
% figure;
% plot(f, angle(Srw1)*180/pi); grid;
% set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
% xlabel('f (cycles/sample)');
% ylabel('Phase (degrees)');
% title('Power spectral estimates Srw1(f)');