%% plot VCO Phase Noise

figure();
semilogx(f, 10*log10(VCO_n1));
hold on; grid on;
semilogx(f, 10*log10(VCO_n2));
semilogx(f, 10*log10(VCO_n3));
semilogx(f, 10*log10(Sn_VCO_OL), LineWidth=3);

legend('Flicker', 'Thermal', 'Floor', 'Total');
xlabel('Offset Frequency (Hz)');
ylabel('Phase Noise (dBc/Hz)');

%% Xtal phase noise Plotting
figure;
loglog(x, y, 'bo', 'MarkerFaceColor', 'b'); % Original data points
hold on;
loglog(x_new, y_new, 'r-'); % Interpolated data
xlabel('Frequency (Hz)');
ylabel('Noise Level (dBc/Hz)');
legend('Original data', 'Log-Log Interpolated', 'Location', 'Southwest');
grid on;
title('Log-Log Interpolation of Noise Level vs. Frequency');

%% Closing loop for VCO noise

% forward path over loop gain
H_nVCO                                                  = 1 / (1 + GOL);

[H_nVCO_mag, H_nVCO_phase, w]          = bode(H_nVCO, w);

H_nVCO_mag                                          = squeeze(H_nVCO_mag);


Sn_VCO_CL                                               = Sn_VCO_OL' .* H_nVCO_mag.^2;



