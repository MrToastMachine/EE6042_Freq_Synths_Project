close all;
clear all;

% Definitions:
s = tf('s');

f = logspace(0, 7, 10000);
w = 2*pi*f;

% Ref & VCO Freqs:
fref = 25e6;
fvco = 2e9;

% XOR PD Supply
VDD = 1;

Kpd = VDD/pi;

% Loop Filter:
wp = 2*pi*365.76e3;
Hlf = 1/(1+s/wp);

% VCO:
Kv = 2*pi*63.75e6;

% F/B Divider:
N = fvco/fref;

% Open & Closed Loop TFs:
Gol = Kpd * Hlf * Kv/s * (1/N);

Gcl = Gol / (1 + Gol);

[Gol_mag, Gol_phase, wout] = bode(Gol,w);

Gol_mag = squeeze(Gol_mag);
Gol_phase = squeeze(Gol_phase);

[Gcl_mag, Gcl_phase, wout] = bode(Gcl,w);

Gcl_mag = squeeze(Gcl_mag);
Gcl_phase = squeeze(Gcl_phase);

% Bode plots:
figure(1);
subplot(2,1,1);
semilogx(w/(2*pi), 20*log10(Gol_mag), 'LineWidth', 2, 'DisplayName', 'Mag Gol'); grid on; hold on;
title('Gol Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');


subplot(2,1,2);
semilogx(w/(2*pi), Gol_phase, 'LineWidth', 2, 'DisplayName', 'Phase Gol'); grid on; hold on;
title('Gol Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');

figure(1);
subplot(2,1,1);
semilogx(w/(2*pi), 20*log10(Gcl_mag), 'LineWidth', 2, 'DisplayName', 'Mag Gcl'); grid on; hold on;
title('Gcl Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Show');

subplot(2,1,2);
semilogx(w/(2*pi), Gcl_phase, 'LineWidth', 2, 'DisplayName', 'Phase Gol'); grid on; hold on;
title('Gcl Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Frequency (Hz)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');
legend('Show');