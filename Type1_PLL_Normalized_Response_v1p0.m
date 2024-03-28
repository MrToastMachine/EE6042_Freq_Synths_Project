close all;
clear all;

%% Definitions
% Define frequency space:
f = logspace(-6, 2, 10000);
w = 2*pi*f;

% Unity gain freq for 'no-filter' case:
w0 = 1;

% Sweep range for damping factor:
zeta = 0.1:0.1:2;

% Laplace varaiable:
s = tf('s');

%% Main For Loop:
% Loop over zeta
% For each zeta, calc:
% * wp
% * wn
% * Gol(s)
% * Gcl(s)
% * PM
% * -3dB CL BW
%
% Plot:
% * Bode plots for Gol, Gcl
% * Closed loop step response
% * Find poles of the system, and plot root locus

for ii = 1:length(zeta)
    wp(ii) = zeta(ii)^2 * 4 * w0;
    wn(ii) = sqrt(wp(ii) * w0);
    
    Gol = w0/s*wp(ii)/(s+wp(ii));
    
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
    semilogx(w, 20*log10(Gol_mag), 'LineWidth', 2); grid on; hold on;
    title('Gol Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Zeta = 0.1','Zeta = 0.2','Zeta = 0.3','Zeta = 0.4','Zeta = 0.5','Zeta = 0.6','Zeta = 0.7','Zeta = 0.8','Zeta = 0.9','Zeta = 1.0','Zeta = 1.1','Zeta = 1.2','Zeta = 1.3','Zeta = 1.4','Zeta = 1.5','Zeta = 1.6','Zeta = 1.7','Zeta = 1.8','Zeta = 1.9','Zeta = 2.0');
    
    subplot(2,1,2);
    semilogx(w, Gol_phase, 'LineWidth', 2); grid on; hold on;
    title('Gol Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Zeta = 0.1','Zeta = 0.2','Zeta = 0.3','Zeta = 0.4','Zeta = 0.5','Zeta = 0.6','Zeta = 0.7','Zeta = 0.8','Zeta = 0.9','Zeta = 1.0','Zeta = 1.1','Zeta = 1.2','Zeta = 1.3','Zeta = 1.4','Zeta = 1.5','Zeta = 1.6','Zeta = 1.7','Zeta = 1.8','Zeta = 1.9','Zeta = 2.0');
    
    figure(2);
    subplot(2,1,1);
    semilogx(w, 20*log10(Gcl_mag), 'LineWidth', 2); grid on; hold on;
    title('Gcl Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Zeta = 0.1','Zeta = 0.2','Zeta = 0.3','Zeta = 0.4','Zeta = 0.5','Zeta = 0.6','Zeta = 0.7','Zeta = 0.8','Zeta = 0.9','Zeta = 1.0','Zeta = 1.1','Zeta = 1.2','Zeta = 1.3','Zeta = 1.4','Zeta = 1.5','Zeta = 1.6','Zeta = 1.7','Zeta = 1.8','Zeta = 1.9','Zeta = 2.0');
    
    subplot(2,1,2);
    semilogx(w, Gcl_phase, 'LineWidth', 2); grid on; hold on;
    title('Gcl Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Zeta = 0.1','Zeta = 0.2','Zeta = 0.3','Zeta = 0.4','Zeta = 0.5','Zeta = 0.6','Zeta = 0.7','Zeta = 0.8','Zeta = 0.9','Zeta = 1.0','Zeta = 1.1','Zeta = 1.2','Zeta = 1.3','Zeta = 1.4','Zeta = 1.5','Zeta = 1.6','Zeta = 1.7','Zeta = 1.8','Zeta = 1.9','Zeta = 2.0');
    
    % Step Response Plot
    t=0:1e-3:100;

    phase_out = step(Gcl,t); %%step response of the system.
    
    figure(3);
    plot(t,phase_out, 'LineWidth', 2); grid on; hold on
    title ('Step Response of PLL for a Unit Step Input Phase Change', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (rad)', 'FontSize', 14, 'FontWeight', 'bold');
    legend('Zeta = 0.1','Zeta = 0.2','Zeta = 0.3','Zeta = 0.4','Zeta = 0.5','Zeta = 0.6','Zeta = 0.7','Zeta = 0.8','Zeta = 0.9','Zeta = 1.0','Zeta = 1.1','Zeta = 1.2','Zeta = 1.3','Zeta = 1.4','Zeta = 1.5','Zeta = 1.6','Zeta = 1.7','Zeta = 1.8','Zeta = 1.9','Zeta = 2.0');
    
    P = pole(minreal(Gcl));
    figure(5);
    plot(real(P(1)),imag(P(1)), 'x'); hold on, grid on;
    plot(real(P(2)),imag(P(2)), 'x'); hold on, grid on;
    
    Aol_dB = 20*log10(Gol_mag);
    index1 = find(Aol_dB <= 0);
    fu_array(ii) = w(index1(1))/(2*pi);
    
    % Phase Margin
    PM_array(ii) = 180 + (unwrap(Gol_phase(index1(1))));
    
    % Find the -3dB CL BW with Ideal LF:
    Acl_dB = 20*log10(Gcl_mag);
    index2 =find(Acl_dB <= -3);
    fcl_bw_array(ii) = w(index2(1))/(2*pi);
end

%% Plot PM vs Damping Factor:
figure(44);
yyaxis left;
plot(zeta, PM_array, 'LineWidth', 2); grid on;
xlabel('Damping Factor (zeta)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Phase Margin (deg)', 'FontSize', 14, 'FontWeight', 'bold');

yyaxis right;
plot(zeta, fcl_bw_array*2*pi, 'LineWidth', 2); grid on;
ylabel('w_-_3_d_B/w_0', 'FontSize', 14, 'FontWeight', 'bold');

PM = 180/pi * atan( 2*zeta ./ sqrt(-2*zeta.^2 + sqrt(1+4*zeta.^4)) );

figure(4);
plot(zeta, PM, 'LineWidth', 2); grid on;
title('Phase Margin vs Damping Factor', 'FontSize', 14, 'FontWeight', 'bold');
xlabel('Damping Factor (zeta)', 'FontSize', 14, 'FontWeight', 'bold');
ylabel('Phase Margin (deg)', 'FontSize', 14, 'FontWeight', 'bold');

%% Critically Damped Case:
zeta = 1/sqrt(2);
    wp = zeta^2 * 4*w0;
    wn = sqrt(wp * w0);
    
    Gol = w0/s*wp/(s+wp);
    
    Gcl = Gol / (1 + Gol);
    
    [Gol_mag, Gol_phase, wout] = bode(Gol,w);

    Gol_mag = squeeze(Gol_mag);
    Gol_phase = squeeze(Gol_phase);

    [Gcl_mag, Gcl_phase, wout] = bode(Gcl,w);

    Gcl_mag = squeeze(Gcl_mag);
    Gcl_phase = squeeze(Gcl_phase);

    figure(10);
    subplot(2,1,1);
    semilogx(w, 20*log10(Gcl_mag), 'LineWidth', 2); grid on; hold on;
    title('Gcl Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');
    
    subplot(2,1,2);
    semilogx(w, Gcl_phase, 'LineWidth', 2); grid on; hold on;
    title('Gcl Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');
    
    figure(11);
    subplot(2,1,1);
    semilogx(w, 20*log10(Gol_mag), 'LineWidth', 2); grid on; hold on;
    title('Gol Magnitude Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Magnitude Response (dB)', 'FontSize', 14, 'FontWeight', 'bold');
        
    subplot(2,1,2);
    semilogx(w, Gol_phase, 'LineWidth', 2); grid on; hold on;
    title('Gol Phase Response', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('w/w_0 Normalized Frequency', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (deg)', 'FontSize', 14, 'FontWeight', 'bold');

    % Step response of the system
    phase_out = step(Gcl,t); 
    
    figure(30);
    plot(t,phase_out, 'LineWidth', 2); grid on; hold on
    title ('Step Response of PLL for a Unit Step Input Phase Change', 'FontSize', 14, 'FontWeight', 'bold');
    xlabel('Time (s)', 'FontSize', 14, 'FontWeight', 'bold');
    ylabel('Phase Response (rad)', 'FontSize', 14, 'FontWeight', 'bold');