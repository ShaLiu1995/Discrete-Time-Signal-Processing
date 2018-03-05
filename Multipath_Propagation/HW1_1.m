h = [1, 0, 0, 0, -0.5];
n = 0 : 4;
stem(n, h);
xlim([-1, 5]);
ylim([-1, 1.5]);
axis([-1 5 -1 1.5]);
xlabel('n');
ylabel('h[n]');
title('h[n]');
grid;

syms z
H_z = zTransform(h);
[num, den] = numden(H_z);
p = double(solve(den == 0, z));
z = double(solve(num == 0, z));
figure;
zplane(p, z);
title('Zero-pole plot of H(z)');
grid;

NFFT = 256;
h_256 = [h, zeros(1, NFFT - 5)];


fs = 1;
f = fs * (0 : NFFT - 1) / NFFT;
H_256 = fft(h_256, NFFT);
magH = abs(H_256);
phaH = angle(H_256);
figure();
plot(f, magH);
title('|H[k]| vs f');
xlabel('Normalized Frequency')       
ylabel('Amplitude');
axis([0 0.5*fs 0 2.5]);
grid;

figure;
plot(f, phaH);
title('arg(H[k]) vs f');
xlabel('Normalized Frequency')       
ylabel('Phase/Radians');
axis([0 0.5*fs -pi pi]);
grid;
