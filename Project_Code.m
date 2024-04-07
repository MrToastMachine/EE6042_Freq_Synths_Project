%% Loop Design

s = tf('s');
w = 2*pi*logspace(1,8,1000);
f = w/(2 * pi);

f_xtal = 38.4e6;
f_vco = 1.8433e9

N = f_vco / f_xtal;


I_cp = 4e-3 % 4mA
I_cp0 = 4e-3;

K_pd = I_cp / (2*pi);
K_v = 2 * pi * 30e6;

R1 = 1000;
C1 = 1e-6;

H_lf = (1+s*R1*C1)/(s*C1);


GOL = K_pd * H_lf * (K_v/s) * (1/N);
GCL = GOL / (1 + GOL);

bode(GCL)

%% VCO Phase Noise profile

fc_flicker_Hz = 1e3; % flicker noise at 1kHz
fc_thermal_Hz = 2e6; % thermal noise at 2MHz

PN_flicker_dBcHz = -44.5;
PN_thermal_dBcHz = -130;

VCO_floor_dBcHz = -150;

VCO_flicker_spot_freq         = fc_flicker_Hz;
VCO_flicker_dBcHz              = PN_flicker_dBcHz;

VCO_thermal_spot_freq       = fc_thermal_Hz;
VCO_thermal_dBcHz            = PN_thermal_dBcHz;

VCO_n1 = 10.^((VCO_flicker_dBcHz - 30.*log10(f) + 30.*log10(VCO_flicker_spot_freq))./10);
VCO_n2 = 10.^((VCO_thermal_dBcHz - 20.*log10(f) + 20.*log10(VCO_thermal_spot_freq))./10);
VCO_n3 = 10^(VCO_floor_dBcHz/10) * ones(size(f));

Sn_VCO_OL = (VCO_n1 + VCO_n2 + VCO_n3);


%% Closing loop for VCO noise

H_nVCO                                                  = 1 / (1 + GOL);

[H_nVCO_mag, H_nVCO_phase, w]          = bode(H_nVCO, w);

H_nVCO_mag                                          = squeeze(H_nVCO_mag);


Sn_VCO_CL                                               = Sn_VCO_OL' .* H_nVCO_mag.^2;

%% plot VCO Phase Noise OL and CL

figure();
semilogx(f, 10*log10(Sn_VCO_OL), LineWidth=2);
hold on; grid on;
semilogx(f, 10*log10(Sn_VCO_CL), LineWidth=2);

legend('Open Loop', 'Closed Loop');
xlabel('Offset Frequency (Hz)');
ylabel('Phase Noise (dBc/Hz)');

%% Xtal Phase Noise Interpolation

% known xtal phase noise points
x = [100, 1e3, 10e3, 100e3, 1e6, 10e6]; % Frequencies in Hz
y = [-110, -140, -150, -150, -150, -150]; % Noise levels in dBc/Hz

% New frequencies where you want to interpolate noise levels
x_new = logspace(log10(min(x)), log10(max(x)), 1000); % 1000 new points between min and max of x

% Convert data to log scale for both x and y
log_x = log10(x);
log_y = log10(y);

% Convert the new x values to log scale
log_x_new = log10(x_new);

% Perform linear interpolation in the log-log space
log_y_new = interp1(log_x, log_y, log_x_new, 'pchip');

% Convert interpolated values back to the original scale
temp = 10.^log_y_new;
Sn_Xtal_OL = 10.^temp;

% forward path over loop gain
H_nXtal = (K_pd*H_lf*(K_v/s))/(1 + GOL);

[H_nXtal_mag, H_nXtal_phase, w] = bode(H_nXtal, w);

H_nXtal_mag = squeeze(H_nXtal_mag);

Sn_Xtal_CL = Sn_Xtal_OL' .* H_nXtal_mag.^2;

semilogx(f, log10(Sn_Xtal_OL));


%% CP Noise

i_2n = (I_cp / I_cp0) * (2.5e-21);

% forward path over loop gain
CPn_noise = (H_lf * (K_v / s)) / (1 + GOL);

[CP_noise_mag, CP_noise_phase, w] = bode(CPn_noise);

CP_noise_mag = squeeze(CP_noise_mag);

Sn_CP_CL = CP_noise_mag^2 * i_2n * 0.5;

bode(Sn_CP_CL);


%% Noise Mask
close all
mask_freqs = cat(2, linspace(1e4, 2e6 -1, 100), linspace(2e6, 1e7,100));
mask_vals = cat(2,-100*ones(100), -120*ones(100));

figure()
semilogx(f, 10*log10(Sn_VCO_CL), 'LineWidth', 2);
grid on; hold on;
semilogx(mask_freqs, mask_vals, '--', 'LineWidth', 2)