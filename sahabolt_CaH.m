function [ sahabolt ] = sahabolt_CaH(tempCaH,k,p_CaH,ionCaH,lvl)
% compute Saha-Boltzmann populaton n_(r,s)/N for any level r,s of E,
% input: temperature, boltzmann's constant, electron pressure, ion stage, ionization level

saha = saha_CaH(tempCaH,p_CaH);
sahabolt = saha(1,ionCaH) * boltz_CaH(tempCaH,k,ionCaH,lvl); 
    %using a designated vector element of saha according to ion stage

end
