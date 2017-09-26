function [ n_rs ] = boltz_CaH(tempCaH,k,r,s)
% compute Boltzmann population for level r,s of Schadee element E
% input: temperature, ion stage r, level s
% output: relative level population n_(r,s)/N_r
    
    u = partfunc_CaH(tempCaH,k);     % get partition functions
    n_rs = 1 / u(r) * exp (-(s-1) / (k*tempCaH)); % Boltzmann distribution
    
end

