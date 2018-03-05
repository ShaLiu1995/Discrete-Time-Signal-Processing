% 1
N = 1024;
n = 0:1023;
A = sqrt(2);
signal = A * sin(2*pi*0.125*n);
noise = randn(1,1024);
x = signal + noise;
figure;
plot(n, x);
grid;
xlabel('x(n)');
ylabel('n');
title ('Time series of the sinusoid signal plus the Guassain white noise');
xlim([0,1024]);

% 2
fs = 1;
for NFFT = [128, 256, 512, 1024]
    [pxx,f] = pwelch(x(1:NFFT), NFFT, 0, NFFT, fs, 'centered');
    figure; grid;
    plot(f, 10*log10(pxx));
    set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
    title(['Power Sprectral Estimatioin , N = NFFT = ' num2str(NFFT) ]);
    ylim([-30 30]);
    xlabel('f (cycles/sample)');
    ylabel("Magnitude(dB)");
    grid;
    
    figure; grid;
    plot(f, 10*log10(pxx));
    set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
    title(['Power Sprectral Estimatioin (Zoom in) , N = NFFT = ' num2str(NFFT) ]);
    ylim([-5 5]);
    xlabel('f (cycles/sample)');
    ylabel("Magnitude(dB)");
    grid;
end

3
M = 128;
for overlap = [0, 64, 96]
    window = hamming(M);
    [pxx,f] = pwelch(x, M, overlap, M, fs,'centered');
    U = sum(window.^2) / M;
    figure; grid;
    plot(f, 10*log10(pxx));
    ylim([-5 20]);
    set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
    xlabel('f (cycles/sample)');
    ylabel('Magnitude (dB)');
    title(['Power Spectral Estimation, N=1024, M=128, Overlap=' num2str(overlap)]);
    grid;
    
    figure; grid;
    plot(f, 10*log10(pxx));
    set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
    ylim([-5 5]);
    xlabel('f (cycles/sample)');
    ylabel('Magnitude (dB)');
    title(['Power Spectral Estimation (Zoom in), N=1024, M=128, Overlap=' num2str(overlap)]);
    grid;
end

% 4
% segment the time series
x_reshape = reshape(x,[128, 8])';
% calculate the coherent average
x_avg = sum(x_reshape, 1)./8;
n = 0:127;
figure;
plot(n, x_avg);
grid;
xlabel('n');
ylabel('value');
title ('Coherently averaged time series');
xlim([0 128]);
ylim([-2 2]);

NFFT = 128;
M = 128;
window = hamming(M);
U = sum(window.^2) / M;
f = -0.5 : 1/NFFT : 0.5-1/NFFT;
X = fft(x_avg.*window' ,NFFT);
X = fftshift(X);
figure;
plot(f, 10*log10(abs(X).^2/(fs*M*U)));
grid;
% ylim([-5 20]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power Spectral Estimation (coherently averaging before FFT)');

X2_matrix = zeros(8, 128);
for i=1:8
   X2_matrix(i,:) = fftshift(fft(x_reshape(i,:).*window', NFFT)); 
end
X2 = sum(X2_matrix, 1)./8;
figure;
plot(f, 10*log10(abs(X2).^2/(fs*M*U)));
grid;
% ylim([-5 20]);
set(gca,'xtick', [-0.5 -0.25 0 0.25 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('Power Spectral Estimation (FFT before coherently averaging)');



