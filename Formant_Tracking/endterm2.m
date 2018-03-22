load -ascii jcnwwa_single_col.txt

wwa = jcnwwa_single_col';
NFFT = 256;
beta = 7.85;
Fs = 10000;
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
window = kaiser(NFFT, beta);
MU = sum(window.^2);
p = 14;

% A & B
n = 1:15360;
t = n/Fs;
figure;
plot(t, wwa); grid;
title('Fig. 1. Waveform of entire speech phrase');
xlabel('Time (secs)'); 
ylabel('Amplitude');
xlim([0 1.536]);
% figure;
% spectrogram(wwa/max(wwa),256,128,NFFT,Fs,'yaxis');
% title('Fig. 2. Spectrogram of entire speech phrase');

n = 1:10240;
t = n/Fs;
wwa_extract = wwa(4400:14639);
figure;
plot(t, wwa_extract); grid;
title('Fig. 2. Waveform of extracted speech phrase');
xlabel('Time (secs)'); 
ylabel('Amplitude');
xlim([0 1.024]);
% figure;
% spectrogram(wwa/max(wwa),256,128,NFFT,Fs,'yaxis');
% title('Fig. 4. Spectrogram of extracted speech phrase');

% C
wwa_segment = zeros(256, 79);
for k = 1:79
    wwa_segment(:,k) = wwa_extract((k-1)*128+1 : (k+1)*128);
end

WWA_FFT = zeros(256, 79);
for k = 1:79
     Xk = fftshift(fft(wwa_segment(:,k)'.*window',NFFT));
     WWA_FFT(:,k) = 10*log10((abs(Xk).^2)/(Fs*MU))';
end
C = WWA_FFT(129:256,:);
x = (1:79)*128/Fs;
y = (0:0.5*Fs)/1000;
figure;
imagesc(x,y,C), axis xy;
title('Fig. 3. Conventional Spectrogram of extracted speech phrase');
xlabel("Time (secs)");
ylabel("Frequency (kHz)");
cb = colorbar;
ylabel(cb, 'Magniture (dB)');


WWA_LPC = zeros(256, 79);
for k = 1:79
    a_hat = lpc(wwa_segment(:,k)'.*window', p);
    A_HAT = fftshift(fft(a_hat, NFFT));
    WWA_LPC(:,k) = 10*log10(1./(abs(A_HAT).^2))';
end
C = WWA_LPC(129:256,:);
x = (1:79)*128/Fs;
y = (0:0.5*Fs)/1000;
figure;
imagesc(x,y,C), axis xy;
title('Fig. 4. LPC Spectrogram of extracted speech phrase');
xlabel("Time (secs)");
ylabel("Frequency (kHz)");
cb = colorbar;
ylabel(cb, 'Magniture (dB)');

% d
figure; 
for k = 1:8
    startIndex = 5000+(k-1)*64;
    xn = wwa_extract(startIndex:startIndex+255);
    Xk = fftshift(fft(xn.*window',NFFT));
    subplot(1,2,1);
    plot(f/1000, 10*log10((abs(Xk).^2)/(Fs*MU)) - 50*(k-1));
    xlabel('Frequency (kHz)');
    title('Conventional FFT analysis');
    hold on;
    plot(f/1000, -50*(k-1)*ones(1, 256),'k--');
    hold on;
end
grid on;
% pbaspect([1 2 1]);
xlim([0 0.5]*Fs/1000);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*Fs/1000);
set(gca,'ytick', []);

for k = 1:8
    startIndex = 5000+(k-1)*64;
    xn = wwa_extract(startIndex:startIndex+255);
    a_hat = lpc(xn.*window', p);
    A_HAT = fftshift(fft(a_hat, NFFT));
    subplot(1,2,2);
    plot(f/1000, 10*log10(1./(abs(A_HAT).^2)) - 50*(k-1));
    xlabel('Frequency (kHz)');
    title('LPC analysis');
    hold on;
    plot(f/1000, -50*(k-1)*ones(1, 256),'k--');
    hold on;
end
grid on;
% pbaspect([1 2 1]);
xlim([0 0.5]*Fs/1000);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*Fs/1000);   
set(gca,'ytick', []);
suptitle('Fig. 5. Waterfall plot pair (stationary)');

figure; 
for k = 1:8
    startIndex = 2500+(k-1)*64;
    xn = wwa_extract(startIndex:startIndex+255);
    Xk = fftshift(fft(xn.*window',NFFT));
    subplot(1,2,1);
    plot(f/1000, 10*log10((abs(Xk).^2)/(Fs*MU)) - 50*(k-1));
    xlabel('Frequency (kHz)');
    title('Conventional FFT analysis');
    hold on;
    plot(f/1000, -50*(k-1)*ones(1, 256),'k--');
    hold on;
end
grid on;
% pbaspect([1 2 1]);
xlim([0 0.5]*Fs/1000);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*Fs/1000);
set(gca,'ytick', []);

for k = 1:8
    startIndex = 2500+(k-1)*64;
    xn = wwa_extract(startIndex:startIndex+255);
    a_hat = lpc(xn.*window', p);
    A_HAT = fftshift(fft(a_hat, NFFT));
    subplot(1,2,2);
    plot(f/1000, 10*log10(1./(abs(A_HAT).^2)) - 50*(k-1));
    xlabel('Frequency (kHz)');
    title('LPC analysis');
    hold on;
    plot(f/1000, -50*(k-1)*ones(1, 256),'k--');
    hold on;
end
grid on;
% pbaspect([1 2 1]);
xlim([0 0.5]*Fs/1000);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]*Fs/1000);   
set(gca,'ytick', []);
suptitle('Fig. 6. Waterfall plot pair (rapidly-evolving)');

% e
a = 1;
b = [1 -1];
wwa_filter = filter(b,a,wwa);
wwa_segment_filter = zeros(256, 79);
for k = 1:119
    wwa_segment_filter(:,k) = wwa_filter((k-1)*128+1 : (k+1)*128);
end

WWA_LPC_filter = zeros(256, 119);
for k = 1:119
    a_hat = lpc(wwa_segment_filter(:,k)'.*window', p);
    A_HAT = fftshift(fft(a_hat, NFFT));
    WWA_LPC_filter(:,k) = 10*log10(1./(abs(A_HAT).^2))';
end

wwa_lpc_peak_locs = zeros(3, 119);
for k = 1:119
    [peaks,locs] = findpeaks(WWA_LPC_filter(129:256,k));
    wwa_lpc_peak_locs(:,k) = locs(1:3)*Fs/NFFT/1000;
end

t = (1:119)*128/Fs;
figure;
plot(t, wwa_lpc_peak_locs(1,:),'r+');
hold on;
plot(t, wwa_lpc_peak_locs(2,:),'b+');
hold on;
plot(t, wwa_lpc_peak_locs(3,:),'k+');
hold on;
grid on;
legend('F_1','F_2','F_3');
title('Fig. 7. Peak-picked plot of formant trajectories');
xlabel('Time (secs)'); 
ylabel('Frequency (kHz)');
xlim([0.44 1.536]);
    