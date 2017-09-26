function [ nstage_rel ] = saha_CaH(tempCaH,p_CaH)
% compute Saha population fraction N_r/N for Schadee element E
% input: temperature, electron pressure [dyne/cm^2], ion stage = r

% Constants
k = 8.61734e-5;             % Boltzmann's constant [eV/K]
kerg = 1.380658e-16;        % Boltzmann's constant [erg/K]
h = 6.62607e-27;            % Planck's constant [erg*s]
M_e = 9.109390e-28;         % electron mass [g]
N_e = p_CaH ./ (kerg*tempCaH);   % electron density [cm^-3] (1 erg = 1 dyne*cm)

chiCaH = [6.113 11.871 50.91 67.15];  % Ca ionization energies [eV]

u = partfunc_CaH(tempCaH,k);     % get partition functions U[1]...U[4]
u(end+1) = 2;               % adding a new element to the array, estimated fifth value to get N_4 as well

sahaconst = ((2*pi*M_e*kerg*tempCaH)/(h^2))^(3/2) * 2/N_e; % the constant terms of the equation

nstage = zeros([1 5]);      % 5 zero-element array
nstage(:,1) = 1;            % setting the first element = 1


for r = 1:4
    nstage(r+1) = nstage(r)*sahaconst * (u(r+1)/u(r)) * exp(-chiCaH(1,r)/(k*tempCaH));
    ntotal = sum(nstage);           % sum all stages = element density
    nstage_rel = nstage ./ ntotal;  % fractions of element density
end


end

