% Example data points
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
y_new = 10.^log_y_new;

% Plotting for visualization
figure;
loglog(x, y, 'bo', 'MarkerFaceColor', 'b'); % Original data points
hold on;
loglog(x_new, y_new, 'r-'); % Interpolated data
xlabel('Frequency (Hz)');
ylabel('Noise Level (dBc/Hz)');
legend('Original data', 'Log-Log Interpolated', 'Location', 'Southwest');
grid on;
title('Log-Log Interpolation of Noise Level vs. Frequency');


