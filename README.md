# ISSA_RADAR_SYSTEM

clc;
clear;
close all;

%% -------------------- CONSTANTS --------------------
c = 3e8;                  % Speed of light (m/s)
k = 1.38e-23;             % Boltzmann constant

%% -------------------- RADAR PARAMETERS --------------------
Pt = 0.25e6;              % Peak transmitted power (W)
f = 2.25e9;               % Frequency (Hz)
lambda = c/f;             % Wavelength (m)

G_dB = 30;                % Antenna gain (dB)
G = 10^(G_dB/10);         % Linear gain

tau = 8e-6;               % Pulse width (s)
B = 1e6;                  % Bandwidth (Hz)
NF_dB = 6;                % Noise figure (dB)
NF = 10^(NF_dB/10);       

T = 300;                  % Noise temperature (K)
L_dB = 5;                 % System losses (dB)
L = 10^(L_dB/10);

sigma = 1.5;              % Target RCS (m^2)
Pfa = 1e-6;               % Probability of false alarm

%% -------------------- GEOMETRY --------------------
ht = 1500;                % Target height (m)
hr = 30;                  % Radar height (m)
R = linspace(5e3,200e3,500); % Range (m)

%% -------------------- MULTIPATH (2-RAY MODEL) --------------------
d1 = sqrt(R.^2 + (ht-hr).^2);
d2 = sqrt(R.^2 + (ht+hr).^2);
phi = 2*pi*(d2-d1)/lambda;
Gamma = -1;               % Perfectly reflecting flat earth

F = abs(1 + Gamma.*exp(-1j*phi)); % Pattern propagation factor

%% -------------------- RADAR RANGE EQUATION --------------------
Pr = (Pt * G^2 * lambda^2 * sigma .* F.^2) ./ ...
     ((4*pi)^3 * R.^4 * L);

%% -------------------- NOISE POWER --------------------
Pn = k * T * B * NF;

SNR = Pr ./ Pn;
SNR_dB = 10*log10(SNR);

%% -------------------- (a) SNR vs Range --------------------
figure;
plot(R/1e3, SNR_dB, 'LineWidth', 2);
xlabel('Range (km)');
ylabel('SNR (dB)');
title('SNR vs Range (Target height = 1500 m)');
grid on;

%% -------------------- (b) Probability of Detection --------------------
alpha = sqrt(-2*log(Pfa));
Pd = marcumq(sqrt(2*SNR), alpha);

figure;
plot(R/1e3, Pd, 'LineWidth', 2);
xlabel('Range (km)');
ylabel('Probability of Detection');
title('Pd vs Range (Target height = 1500 m)');
grid on;

%% -------------------- (c) Azimuth Coverage --------------------
az = -60:0.2:60;          % Azimuth sweep (deg)
G_az = G * abs(cosd(az));
Pr_az = (Pt * G_az.^2 * lambda^2 * sigma) ./ ...
        ((4*pi)^3 * (80e3)^4 * L);
SNR_az = Pr_az / Pn;
Pd_az = marcumq(sqrt(2*SNR_az), alpha);

figure;
plot(az, Pd_az, 'LineWidth', 2);
yline(0.75,'r--');
xlabel('Azimuth angle (deg)');
ylabel('Probability of Detection');
title('Azimuth Coverage (Pd = 0.75)');
grid on;

%% -------------------- (d) Elevation Coverage --------------------
el = 0:0.1:30;            % Elevation sweep (deg)
G_el = G * abs(cosd(el));
Pr_el = (Pt * G_el.^2 * lambda^2 * sigma) ./ ...
        ((4*pi)^3 * (80e3)^4 * L);
SNR_el = Pr_el / Pn;
Pd_el = marcumq(sqrt(2*SNR_el), alpha);

figure;
plot(el, Pd_el, 'LineWidth', 2);
yline(0.75,'r--');
xlabel('Elevation angle (deg)');
ylabel('Probability of Detection');
title('Vertical Coverage (Pd = 0.75)');
grid on;

%% -------------------- (e) SNR vs Range for different Pt --------------------
Pt_values = [0.1 0.25 0.5]*1e6;

figure; hold on;
for i = 1:length(Pt_values)
    Pr_temp = (Pt_values(i)*G^2*lambda^2*sigma)./((4*pi)^3*R.^4*L);
    SNR_temp = Pr_temp/Pn;
    plot(R/1e3,10*log10(SNR_temp),'LineWidth',2);
end
xlabel('Range (km)');
ylabel('SNR (dB)');
legend('0.1 MW','0.25 MW','0.5 MW');
title('SNR vs Range for different Peak Power');
grid on;

%% -------------------- (f) SNR vs Range for different RCS --------------------
sigma_values = [0.5 1.5 5];

figure; hold on;
for i = 1:length(sigma_values)
    Pr_temp = (Pt*G^2*lambda^2*sigma_values(i))./((4*pi)^3*R.^4*L);
    SNR_temp = Pr_temp/Pn;
    plot(R/1e3,10*log10(SNR_temp),'LineWidth',2);
end
xlabel('Range (km)');
ylabel('SNR (dB)');
legend('0.5 m^2','1.5 m^2','5 m^2');
title('SNR vs Range for different RCS');
grid on;
