% B
N = 16;
h = [1, 0, 0, 0, -0.5];

% b.1
H = fft(h, N);
H1 = 1 ./ H;
fs = 1;
f = fs * (0 : N - 1) / N;


stem(f, abs(H));
title('|H[k]| vs f, N = 16');
xlabel('Normalized Frequency')       
ylabel('|H[k]|');
axis([0 0.5*fs 0 2.5]);
grid;

figure;
stem(f, angle(H));
title('arg(H[k]) vs f, N = 16');
xlabel('Normalized Frequency')       
ylabel('arg(H[k])/Radians');
axis([0 0.5*fs -pi pi]);
grid;

figure();
stem(f, abs(H1));
title('|H1[k]| vs f, N = 16');
xlabel('Normalized Frequency')       
ylabel('|H1[k]|');
axis([0 0.5*fs 0 2.5]);
grid;

figure;
stem(f, angle(H1));
title('arg(H1[k]) vs f, N = 16');
xlabel('Normalized Frequency')       
ylabel('arg(H1[k])/Radians');
axis([0 0.5*fs -pi pi]);
grid;

% b.2
h1 = ifft(H1, N);
figure();
n = 0 : N - 1;
stem(n, h1);
title('h1[n]');
xlabel('n');
ylabel('h1[n]');
axis([-1 N -1 1.5]);
grid;

syms z
H2_z = zTransform(h1);
[num, den] = numden(H2_z);
p = double(solve(den == 0, z));
z = double(solve(num == 0, z));
figure;
zplane(p, z);
title('Zero-pole plot of H1(z)');
grid;

% b.3
NFFT = 256;
h1_256 = [h1 zeros(1, NFFT - N)];
H1_256 = fft(h1_256, NFFT);
f = fs * (0 : NFFT - 1) / NFFT;

figure();
plot(f, abs(H1_256));
title('|H1[k]| vs f, N = 16, NFFT = 256');
xlabel('Normalized Frequency')       
ylabel('|H1[k]|');
axis([0 0.5*fs 0 2.5]);
grid;

figure;
plot(f, angle(H1_256));
title('arg(H1[k]) vs f, N = 16, NFFT = 256');
xlabel('Normalized Frequency')       
ylabel('arg(H1[k])/Radians');
axis([0 0.5*fs -pi pi]);
grid;

% b.4
h2 = ifft(fft(h, 2*N).*fft(h1, 2*N));

% b.5
figure()
t = 0 : 2*N - 1;
stem(t, h2);
title('h2[n]');
xlabel('n')       
ylabel('h2[n]');
axis([-1 2*N -1 1.5]);
grid;

syms z
% Truncate the trailing zeros
lastNonZeroIndex = length(h2);
while abs(h2(lastNonZeroIndex)) < 1e-6
    lastNonZeroIndex = lastNonZeroIndex - 1;
end
H2_z = zTransform(h2(1:lastNonZeroIndex));
[num, den] = numden(H2_z);
p = double(solve(den == 0, z));
z = double(solve(num == 0, z));
figure;
zplane(p, z);
title('Zero-pole plot of H2(z)');
grid;

% b.6
NFFT = 256;
h2_256 = [h2 zeros(1, NFFT - N*2)];
H2_256 = fft(h2_256, NFFT);
f = fs * (0 : NFFT - 1) / NFFT;

figure();
plot(f, abs(H2_256));
title('|H2[k]| vs f, 2*N = 32');
xlabel('Normalized Frequency')       
ylabel('|H2[k]|');
xlim([0 0.5]);
ylim([0 2.5]);
axis([0 0.5*fs 0 2.5]);
grid;

figure;
plot(f, angle(H2_256));
title('arg(H2[k]) vs f, 2*N = 32');
xlabel('Normalized Frequency')       
ylabel('arg(H2[k])/Radians');
axis([0 0.5*fs -pi pi]);
grid;

