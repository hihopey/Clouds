clear all
% radar frequency (Hz);
freq = linspace(0.1, 200e9, 1000);
c = 3e11; %speed of light in mm
% radar wavelength (mm);
lambda = c./freq;
% temperature of water in celsius
T_w = [-10 0 10 20];

% Range of normalize diameters (pi D / lambda) to plot
D_n = linspace(0.01, 10, 501);
deltaD = 10;

% Defining variables
nt = length(T_w); % nt = the number of elements in the temperature of the water
nl = length(lambda); % nl = the number of elements in lambda, the radar wavelengths
nd = length(D_n); % nd = the number of elements in the range of normalized diameters

%Fill n_w with values from the interpolated excel equations for each
%temperature
n_w = zeros (nl, nt);
for n = 1: nl
n_w(n, 1) = 9.7979 * exp (-2e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 2) = 9.6233 * exp (-2e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 3) = 9.6649 * exp (-3e-11 .* freq(n));
end
for n = 1: nl
n_w(n, 4) = 8.9104 * exp (-3e-11 .* freq(n));
end

%Fill k_w with values from the interpolated excel equations for each
%temperature
k_w = zeros (nl, nt);
for k = 1: nl
k_w(k, 1) = 0.2666 * reallog((freq(k))) - 3.5211;
end
for k = 1: nl
k_w(k, 2) = -0.3 * reallog((freq(k))) + 9.7807;
end
for k = 1: nl
k_w(k, 3) = 0.7645 * reallog((freq(k))) - 15.309;
end
for k = 1: nl
k_w(k, 4) = 1.0353 * reallog((freq(k))) - 21.881;
end

%Calcualate other components of the complex refractive index needed
m_w = n_w + (i*k_w);
Km = ((m_w.^2)+1)./((m_w.^2)+2);
Km2_w = abs(Km.^2);
ImKm_w = imag(-1 * Km);
 
% set up arrays filling them with zeros
sig_a = zeros(nl, nd); 
sig_s = zeros(nl, nd);
sig_an = zeros(nl, nd);
sig_sn = zeros(nl, nd);
sig_ed = zeros(nl, nd);
att = zeros(nl, 10);
Nd = zeros(nl, nd);
trans = zeros(nl, 10);
LWC = [0.28 0.4 0.3 0.26 0.41 0.28 0.41 0.45 0.3 0.33];
a0 = [5.15 7.73 4.09 2.17 6.04 7 7.74 7.38 3.11 3.51];
u = [10 16.3 2.5 1.6 15.7 7 7.6 4.8 2.8 3.2];
%Ask user which temperature index they want to use
fprintf(1, 'Possible temperature values\n')
fprintf(1, '1: -10C\n')
fprintf(1, '2:  0C\n')
fprintf(1, '3: 10C\n')
fprintf(1, '4: 20C\n')
ti = input('Enter the temperature index (1 - 4): ');
switch ti
 case 1
  ptit = 'Temperature = -10^{\circ}C';
 case 2
  ptit = 'Temperature = 0^{\circ}C';
 case 3
  ptit = 'Temperature = 10^{\circ}C';
 case 4
  ptit = 'Temperature = 20^{\circ}C';
 otherwise
  fprintf('*** Invalid Selection ... exiting!\n');
  return
end
  for i = 1 : 10 
    for lambi = 1 : nl
        for D_n_ind = 1: nd
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            a = (D_n(D_n_ind))/2;
            k = (2 * pi)/(lambda(lambi));
            x = k * a;
            m = m_w(lambi, ti);
            result = Mie(m, x);
            E = result(:, 4) .* (pi*(a^2));
            sig_e = sum(E);
            N = (u(i)^(u(i) + 1))/(gamma(u(i) + 1) * (a0(i)^(u(i) + 1))); 
            Nd = (N * (a^u(i))) * exp(-u(i)*(a/a0(i))); 
            att(lambi, i) = att(lambi, i) + (4.34e3 * Nd * sig_e * deltaD); % calculates attenuation
            trans(lambi, i) = 1 - att(lambi, i); % calculates transmission
        end
    end
  end
    
% Print Attenuation Graph
fntsz = 14;
figure(1)
clf
plot (freq/1e9, att(:, 1), 'b', 'DisplayName', '0.28 g/cm^3')
hold on
plot (freq/1e9, att(:, 2), 'r', 'DisplayName', '0.4 g/cm^3')
plot (freq/1e9, att(:, 3), 'k', 'DisplayName', '0.3 g/cm^3')
plot (freq/1e9, att(:, 4), 'g', 'DisplayName', '0.26 g/cm^3')
plot (freq/1e9, att(:, 5), 'r--', 'DisplayName', '0.41 g/cm^3')
plot (freq/1e9, att(:, 6), 'b', 'DisplayName', '0.28 g/cm^3')
plot (freq/1e9, att(:, 7), 'r--', 'DisplayName', '0.41 g/cm^3')
plot (freq/1e9, att(:, 8), 'k--', 'DisplayName', '0.45 g/cm^3')
plot (freq/1e9, att(:, 9), 'k', 'DisplayName', '0.3 g/cm^3')
plot (freq/1e9, att(:, 10), 'g--', 'DisplayName', '0.33 g/cm^3')
%ylim ([2 10])
lgd = legend('location', 'northwest', 'orientation', 'vertical')
title(lgd, 'LWC')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB / Km)')
title(ptit)
hold off
  
% Print Transmission Graph
figure(2)
clf
plot (freq./1e9, trans(:, 1), 'b', 'DisplayName', '5mm')
hold on
plot (freq/1e9, trans(:, 2), 'r', 'DisplayName', '0.4 g/cm^3')
plot (freq/1e9, trans(:, 3), 'k', 'DisplayName', '0.3 g/cm^3')
plot (freq/1e9, trans(:, 4), 'g', 'DisplayName', '0.26 g/cm^3')
plot (freq/1e9, trans(:, 5), 'r--', 'DisplayName', '0.41 g/cm^3')
plot (freq/1e9, trans(:, 6), 'b', 'DisplayName', '0.28 g/cm^3')
plot (freq/1e9, trans(:, 7), 'r--', 'DisplayName', '0.41 g/cm^3')
plot (freq/1e9, trans(:, 8), 'k--', 'DisplayName', '0.45 g/cm^3')
plot (freq/1e9, trans(:, 9), 'k', 'DisplayName', '0.3 g/cm^3')
plot (freq/1e9, trans(:, 10), 'g--', 'DisplayName', '0.33 g/cm^3')
%ylim ([-10 -2])
lgd = legend('location', 'southwest', 'orientation', 'vertical')
title(lgd, 'LWC')
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel ('Transmission (dB / Km)')
title(ptit)
hold off