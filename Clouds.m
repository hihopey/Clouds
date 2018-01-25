clear all
% radar frequency (Hz);
freq = linspace(0.1, 100e9, 100);
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
att = zeros(nl);
Nd = zeros(nl, nd);
trans = zeros(nl, 6);

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
   
    for lambi = 1 : nl
        for D_n_ind = 1: nd
            D=lambda(lambi)*D_n(D_n_ind)/pi;
            a = D/2;
            k = (2 * pi)/(lambda(lambi));
            x = k * a;
            m = m_w(lambi, ti);
            result = Mie(m, x);
            E = result(:, 4) .* (pi*(a^2));
            sig_e = sum(E);
            Nd = (10^10) * (exp(-100*D)); 
            att(lambi) = att(lambi) + (4.34e3 * Nd * sig_e * deltaD); % calculates attenuation
            trans(lambi) = 1 - att(lambi); % calculates transmission
        end
    end
    
% Print Attenuation Graph
fntsz = 14;
figure(1)
clf
ylim ([2 10])
plot (freq/1e9, att(:, 1), 'b')
hold on
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel('Attenuation (dB / Km)')
title(ptit)
hold off
  
% Print Transmission Graph
fntsz = 14;
figure(2)
clf
ylim ([-10 -2])
plot (freq./1e9, trans(:, 1), 'b')
hold on
set(gca, 'fontsize', fntsz)
xlabel('Frequency (GHz)')
ylabel ('Transmission (dB / Km)')
title(ptit)
hold off