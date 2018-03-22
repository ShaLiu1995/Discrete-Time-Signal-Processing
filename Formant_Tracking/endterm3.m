load -ascii jcnwwa_single_col.txt

wwa = jcnwwa_single_col';
wwa = wwa(4400:14639);
NFFT = 256;
beta = 7.85;
Fs = 10000;
f = -0.5*Fs : Fs/NFFT : (0.5-1/NFFT)*Fs;
window = kaiser(NFFT, beta);
MU = sum(window.^2);
p = 14;

for alpha = -2:0.2:0
    a = 1;
    b = [1 alpha];
    wwa_filter = filter(b,a,wwa);
    wwa_segment_filter = zeros(256, 79);
    for k = 1:79
        wwa_segment_filter(:,k) = wwa_filter((k-1)*128+1 : (k+1)*128);
    end
    WWA_LPC_filter = zeros(256, 79);
    for k = 1:79
        a_hat = lpc(wwa_segment_filter(:,k)'.*window', p);
        A_HAT = fftshift(fft(a_hat, NFFT));
        WWA_LPC_filter(:,k) = 10*log10(1./(abs(A_HAT).^2))';
    end

    wwa_lpc_peak_locs = zeros(3, 79);
    for k = 1:79
        [peaks,locs] = findpeaks(WWA_LPC_filter(129:256,k));
        wwa_lpc_peak_locs(:,k) = locs(1:3)*Fs/NFFT;
    end
    t = (1:79)*128;
    figure;
    plot(t, wwa_lpc_peak_locs(1,:),'r+');
    hold on;
    plot(t, wwa_lpc_peak_locs(2,:),'b+');
    hold on;
    plot(t, wwa_lpc_peak_locs(3,:),'k+');
    hold on;
    legend('F_1','F_2','F_3');
    title('Fig. 7. Peak-picked plot of formant trajectories');
    xlabel('Time (secs)'); 
    ylabel('Frequency (Hz)');
end