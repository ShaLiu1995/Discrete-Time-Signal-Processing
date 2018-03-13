load -ascii bili_single_col.txt
load -ascii bilu_single_col.txt
load -ascii bilng_single_col.txt

aa_i=bili_single_col';
aa_u=bilu_single_col';
aa_ng=bilng_single_col';

aa_i_256 = aa_i(1024:1279);
aa_u_256 = aa_u(1024:1279);
aa_ng_256 = aa_ng(1024:1279);

aa_i_64 = aa_i_256(1:64);
aa_u_64 = aa_u_256(1:64);
aa_ng_64 = aa_ng_256(1:64);


n = 1024:1279;
beta = 7.85;
fs = 10000;
window_256 = kaiser(256, beta);
window_64 = kaiser(64, beta);
MU_256 = sum(window_256.^2);
MU_64 = sum(window_256.^2);
p = 2:2:14;

% A.1
figure;
plot(n, aa_i_256); grid;
title('Figure 1: waveform of /i/');
xlabel('n'); 
ylabel('Amplitude');
xlim([1024 1279]);

figure;
plot(n, aa_u_256); grid;
title('Figure 2: waveform of /u/');
xlabel('n'); 
ylabel('Amplitude');
xlim([1024 1279]);

figure;
plot(n, aa_ng_256); grid;
title('Figure 3: waveform of /ng/');
xlabel('n'); 
ylabel('Amplitude');
xlim([1024 1279]);

% A.2
NFFT = 256;
f = -0.5*fs : fs/NFFT : (0.5-1/NFFT)*fs;

AA_I_256 = fftshift(fft(aa_i_256.*window_256',NFFT));
figure;
plot(f, 10*log10(abs(AA_I_256).^2/(fs*MU_256)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 4: FFT of /i/, 256-point record');

AA_U_256 = fftshift(fft(aa_u_256.*window_256',NFFT));
figure;
plot(f, 10*log10(abs(AA_U_256).^2/(fs*MU_256)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 5: FFT of /u/, 256-point record');

AA_NG_256 = fftshift(fft(aa_ng_256.*window_256',NFFT));
figure;
plot(f, 10*log10(abs(AA_NG_256).^2/(fs*MU_256)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 6: FFT of /ng/, 256-point record');

% A.3
AA_I_64 = fftshift(fft(aa_i_64.*window_64',NFFT));
figure;
plot(f, 10*log10(abs(AA_I_64).^2/(fs*MU_64)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 7: FFT of /i/, 64-point record');

AA_U_64 = fftshift(fft(aa_u_64.*window_64',NFFT));
figure;
plot(f, 10*log10(abs(AA_U_64).^2/(fs*MU_64)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 8: FFT of /u/, 64-point record');

AA_NG_64 = fftshift(fft(aa_ng_64.*window_64',NFFT));
figure;
plot(f, 10*log10(abs(AA_NG_64).^2/(fs*MU_64)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 9: FFT of /ng/, 64-point record');

% B.1 & B.2
% 256-point /i/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_i_256.*window_256', k);
end
figure;
plot(p, Ep); grid;
title('Figure 10: Error variance Ep, /i/, 256-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 16: FFT of 10log(1/|A_{hat}(k)|^2), /i/, 256-point record');

% 256-point /u/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_u_256.*window_256', k);
end
figure;
plot(p, Ep); grid;
title('Figure 11: Error variance Ep, /u/, 256-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2)));
grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 17: FFT of 10log(1/|A_{hat}(k)|^2), /u/, 256-point record');

% 256-point /ng/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_ng_256.*window_256', k);
end
figure;
plot(p, Ep); grid;
title('Figure 12: Error variance Ep, /ng/, 256-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2))); grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 18: FFT of 10log(1/|A_{hat}(k)|^2), /ng/, 256-point record');


% 64-point /i/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_i_64.*window_64', k);
end
figure;
plot(p, Ep); grid;
title('Figure 13: Error variance Ep, /i/, 64-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2))); grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 19: FFT of 10log(1/|A_{hat}(k)|^2), /i/, 64-point record');

% 64-point /u/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_u_64.*window_64', k);
end
figure;
plot(p, Ep); grid;
title('Figure 14: Error variance Ep, /u/, 64-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2))); grid;
xlim([0 0.5*fs]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 20: FFT of 10log(1/|A_{hat}(k)|^2), /u/, 64-point record');

% 64-point /ng/
Ep = zeros(1, 7);
for k = p
    [a_hat, Ep(k/2)] = lpc(aa_ng_64.*window_64', k);
end
figure;
plot(p, Ep); grid;
title('Figure 15: Error variance Ep /ng/, 64-point record');
xlabel('p'); 
ylabel('Ep');
A_hat_k = fftshift(fft(a_hat, NFFT));
figure;
plot(f, 10*log10(1./(abs(A_hat_k).^2))); grid;
xlim([0 0.5]*fs);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*fs);
xlabel('Frequency (Hz)');
ylabel('Magnitude (dB)');
title('Figure 21: FFT of 10log(1/|A_{hat}(k)|^2), /ng/, 64-point record');


