clear all
% radar frequency (GHz);
freq = linspace(100, 200, 1000);
% temperature of water in celsius
T_w = [-10 0 10 20];
% Defining variables
nt = length(T_w); % nt = the number of elements in the temperature of the water
nf = length(freq); % nl = the number of elements in freq, the radar frequencies 
% set up arrays filling them with zeros
att = zeros(nf);
trans = zeros(nf);
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
  t = 263.15;
 case 2
  ptit = 'Temperature = 0^{\circ}C';
  t = 273.15;
 case 3
  ptit = 'Temperature = 10^{\circ}C';
  t = 283.15;
 case 4
  ptit = 'Temperature = 20^{\circ}C';
  t = 293.15;
 otherwise
  fprintf('*** Invalid Selection ... exiting!\n');
  return
end  
  M = input('Liquid water density in the cloud or fog in (g/m^3)?');
  theta = 300 / t;
  fp = 20.2 - 146 * (theta - 1) + 316 * ((theta - 1)^2);
  fs = 39.8 * fp;
  E0 = 77.66 + 103.3 * (theta - 1);
  E1 = 0.0671 * E0;
  E2 = 3.52;
    for f = 1 : nf
            Ep = ((E0-E1)/(1 + ((freq(f)/fp)^2))) + ((E1-E2)/(1 + ((freq(f)/fs)^2))) + E2;
            Ed = ((freq(f) *(E0-E1))/(fp * (1 + ((freq(f)/fp)^2)))) + ((freq(f) *(E1-E2))/(fs * (1 + ((freq(f)/fs)^2))));
            eta = (2 + Ep)/Ed;
            kl = (0.819 * freq(f))/(Ed * (1 + (eta^2)));
            att(f) = kl * M; % calculates attenuation
            trans(f) = 1 - att(f); % calculates transmission
    end