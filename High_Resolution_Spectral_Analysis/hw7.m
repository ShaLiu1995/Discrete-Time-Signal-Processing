% 1.A 1.B
theta_i = [1 2 3 4] * pi/8;
r_i = [0.9896 0.9843 0.9780 0.9686];
z_i = r_i .* exp(1i .* theta_i);
a_n = poly([z_i conj(z_i)]);
% % A_z_recip = 1;
% % syms z_recip
% % for i = 1:4
% %     A_z_recip = A_z_recip * (1-z_i(i)*z_recip) * (1 - conj(z_i(i))*z_recip);
% % end

% a_n = double(coeffs(A_z_recip));
figure;
zplane(1, a_n);
title('Zero-pole plot of H(z)');
grid;

% 1.C
NFFT = 256;
a_n_256 = [a_n zeros(1, NFFT-9)];
f = -0.5 : 1/NFFT : 0.5-1/NFFT;
A_k = fftshift(fft(a_n_256,NFFT));
figure;
plot(f, 10*log10(1./(abs(A_k).^2)));
grid;
xlim([0 0.5]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('FFT of 10log(1/|A(k)|^2)');

% 1.D
n = 0:255;
h_n = filter(1, a_n, [1 zeros(1, 255)]);
figure;
plot(n, h_n); grid;
title('h(n)');
xlabel('n'); 
ylabel('Amplitude');
xlim([0 255]);

H_k = fftshift(fft(h_n, NFFT));
figure();
plot(f, 10*log10(abs(H_k).^2));
grid;
xlim([0 0.5]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('FFT of 10log(|H(k)|^2)');

% 2.A
w_n = randn(1, 1280);
x_n = filter(h_n, 1, w_n);
x_n_256 = x_n(1024:1279);

% 2.B
beta = 7.85;
window = kaiser(NFFT, beta);
fs = 1;
M = 256;
U = sum(window.^2) / M;
X_k = fftshift(fft(x_n_256.*window', NFFT));
figure();
plot(f, 10*log10(abs(X_k).^2/(fs*M*U)));
grid;
xlim([0 0.5]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('FFT of 10log(|X(k)|^2)');

% 2.C
M = 32;
window = kaiser(M, beta);
U = sum(window.^2) / M;

x_n_seg = zeros(15, 256);
index = 1;
for i = 1:15
    x_n_seg(i,:) = [x_n_256(index:index+31).*window' zeros(1, NFFT-32)];
    index = index + 16;
end

X_k_seg = zeros(15, 256);
index = 1;
for i = 1:15
    X_k_seg(i,:) = fftshift(fft(x_n_seg(i,:), NFFT));
end
X_k_avg = sum(abs(X_k_seg).^2)./15;
figure();
plot(f, 10*log10(X_k_avg/(fs*M*U)));
grid;
xlim([0 0.5]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('FFT of 10log(|X(k)|^2 (averaged))');

% 3.A 
window = kaiser(256, beta);
for p = [2, 8, 14]
    [a_hat, g] = lpc(x_n_256.*window', p);
    a_hat_256 = [a_hat zeros(1, NFFT-p-1)];
    A_hat_k = fftshift(fft(a_hat_256, NFFT));
    figure;
    plot(f, 10*log10(1./(abs(A_hat_k).^2)));
    grid;
    xlim([0 0.5]);
    set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
    xlabel('f (cycles/sample)');
    ylabel('Magnitude (dB)');
    title(['FFT of 10log(1/|A_{hat}(k)|^2), p=', num2str(p)]);
    figure;
    zplane(a_hat, 1);
    title(['Zero-pole plot of inverse filter, p=', num2str(p)]);
    grid;  
end

% 3.B
Ep = zeros(1, 7);
for p = 2:2:14
    [a_hat, Ep(p/2)] = lpc(x_n_256.*window', p);
end
p = 2:2:14;
figure;
plot(p, Ep); grid;
title('Error variance Ep');
xlabel('p'); 
ylabel('Ep');


% Repeat II.B
NFFT = 256;
x_n_32 = x_n_256(1:32);
window = kaiser(32, beta);
U = sum(window.^2) / M;
X_k_32 = fftshift(fft(x_n_32.*window', NFFT));
figure();
plot(f, 10*log10(abs(X_k_32).^2/(fs*M*U)));
grid;
xlim([0 0.5]);
set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
xlabel('f (cycles/sample)');
ylabel('Magnitude (dB)');
title('FFT of 10log(|X(k)|^2) (32-point)');

% Repeat III.A
for p = [2, 8, 14]
    [a_hat, g] = lpc(x_n_32.*window', p);
    a_hat_256 = [a_hat zeros(1, NFFT-p-1)];
    A_hat_k = fftshift(fft(a_hat_256, NFFT));
    figure;
    plot(f, 10*log10(1./(abs(A_hat_k).^2)));
    grid;
    xlim([0 0.5]);
    set(gca,'xtick', [0 0.1 0.2 0.3 0.4 0.5]);
    xlabel('f (cycles/sample)');
    ylabel('Magnitude (dB)');
    title(['FFT of 10log(1/|A_{hat}(k)|^2), p=', num2str(p), ' (32-point)']);
    figure;
    zplane(a_hat, 1);
    title(['Zero-pole plot of inverse filter, p=', num2str(p), ' (32-point)']);
    grid;  
end

% Repeat III.B
Ep = zeros(1, 7);
for p = 2:2:14
    [a_hat, Ep(p/2)] = lpc(x_n_32.*window', p);
end
p = 2:2:14;
figure;
plot(p, Ep); grid;
title('Error variance Ep (32-point)');
xlabel('p'); 
ylabel('Ep');
    


