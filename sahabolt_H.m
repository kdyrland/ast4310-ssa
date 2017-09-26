function [ nlevelrel ] = sahabolt_H(temp,p_e,level)
% compute Saha-Boltzmann populaton n_(r,s)/N for any level r,s of E,
% input: temperature, electron pressure, ionization level

k = 8.61734e-5;             % Boltzmann's constant [eV/K]
kerg = 1.380658e-16;        % Boltzmann's constant [erg/K]
kT = k*temp;
kergT = kerg*temp;
h = 6.62607e-27;            % Planck's constant [erg*s]
M_e = 9.109390e-28;         % electron mass [g]
N_e = p_e ./ kergT;   % electron density [cm^-3] (1 erg = 1 dyne*cm)

% Energy levels and weights for hydrogen
nrlvl = 100;                % reasonable partition function cut-off value
g = zeros([2 nrlvl]);       % declarations weights (too many for proton)
chiexc = zeros([2 nrlvl]);  % declaration excitation energies (idem)

for s = 1:nrlvl
    g(1,s) = 2 * s^2;       % statistical weights
    chiexc(1,s) = 13.598*(1-1/(s^2)); % excitation weights
end

g(2,1) = 1;
chiexc(2,1) = 0;


% Partition function
u = zeros([1 2]);
u(1) = 0;

for s = 1:nrlvl
    u(1) = u(1) + g(1,s)*exp(-chiexc(1,s)/kT);
end

u(2) = g(2,1);

% Saha
sahaconst = ((2*pi*M_e*kergT)/(h^2))^(3/2) * 2/N_e;
nstage = zeros([1 2]);
nstage(1) = 1;
nstage(2) = nstage(1) * sahaconst * u(2)/u(1) * exp(-13.598/kT);
ntotal = sum(nstage);           % sum both stages = total hydrogen density

% Boltzmann
nlevel = nstage(1) * g(1,level)/u(1) * exp(-chiexc(1,level)/kT);
nlevelrel = nlevel/ntotal;      % fraction of total hydrogen density

for s = 1:6
    gchi = g(1,s)*exp(-chiexc(1,s)/kT);
%    disp([num2str(s),'  ',num2str(g(1,s)),'  ',num2str(chiexc(1,s)),'  ',num2str(gchi)]) % print
end

%disp(newline)

for s = 1:10:nrlvl      % for s = 1,10,20,...,nrlvl
    gchi = g(1,s)*exp(-chiexc(1,s)/kT);
%    disp([num2str(s),'  ',num2str(g(1,s)),'  ',num2str(chiexc(1,s)),'  ',num2str(gchi)])
end

%disp(newline)

end

