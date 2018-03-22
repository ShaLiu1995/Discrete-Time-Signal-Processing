load -ascii bili_single_col.txt
load -ascii bilu_single_col.txt
load -ascii bilng_single_col.txt
load -ascii jcnwwa_single_col.txt

% a
n = 1:15360;
fs = 10000;
t = n/fs;
wwa = jcnwwa_single_col';
figure;
plot(t, wwa); grid;
title('Figure: waveform');
xlabel('time (sec)'); 
ylabel('Amplitude');
xlim([0 1.536]);

n = 1:10240;
fs = 10000;
t = n/fs;
wwa = wwa(4400:14639);
figure;
plot(t, wwa); grid;
title('Figure: waveform');
xlabel('time (sec)'); 
ylabel('Amplitude');
xlim([0 1.536]);


beta = 7.85;
wwa_e = wwa(0.05*fs : 0.2*fs-1);
window = kaiser(0.15*fs, beta);
MU = sum(window.^2);
WWA_E = fftshift(fft(wwa_e.*window',NFFT));
% figure;
% plot(f, 10*log10(abs(WWA_E).^2/(fs*MU)));
% grid;
% xlim([0 0.5]*fs);
% set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('e');

beta = 7.85;
wwa_ere = wwa(0.4*fs : 0.5*fs-1);
window = kaiser(0.1*fs, beta);
MU = sum(window.^2);
WWA_ERE = fftshift(fft(wwa_ere.*window',NFFT));
% figure;
% plot(f, 10*log10(abs(WWA_ERE).^2/(fs*MU)));
% grid;
% xlim([0 0.5]*fs);
% set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('ere');

beta = 7.85;
wwa_ay = wwa(0.85*fs : 1.0*fs-1);
window = kaiser(0.15*fs, beta);
MU = sum(window.^2);
WWA_AY = fftshift(fft(wwa_ay.*window',NFFT));
% figure;
% plot(f, 10*log10(abs(WWA_AY).^2/(fs*MU)));
% grid;
% xlim([0 0.5]*fs);
% set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('ay');

% b
NFFT = 256;
f = -0.5*fs : fs/NFFT : (0.5-1/NFFT)*fs;
beta = 7.85;
window = kaiser(NFFT, beta);
MU = sum(window.^2);
% spectrogram(wwa/max(wwa),256,128,NFFT,fs,'yaxis');


wwa_segment = zeros(79, 256);
for k = 1:79
    wwa_segment(k,:) = wwa((k-1)*128+1 : (k+1)*128);
end

% WWA_FFT = zeros(79, 256);
% for k = 1:79
%     WWA_FFT(k,:) = fftshift(fft(wwa_segment(k,:).*window',NFFT));
% end


p = 14;
WWA_LPC = zeros(79, 256);
for k = 1:79
    A_hat = lpc(wwa_segment(k,:).*window', p);
    WWA_LPC(k,:) = fftshift(fft(A_hat, NFFT));
end

% [a_hat, ep] = lpc(aa_i_256.*window', p);
% A_hat_k = fftshift(fft(a_hat, NFFT));
% figure;
% plot(f, 10*log10(1./(abs(A_hat_k).^2)));
% grid;
% xlim([0 0.5]*fs);
% set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
% xlabel('Frequency (Hz)');
% ylabel('Magnitude (dB)');
% title('Figure 16: FFT of 10log(1/|A_{hat}(k)|^2), /i/, 256-point record');


