function [ sahabolt ] = sahabolt_E(temp,k,p_e,ion,lvl)
% compute Saha-Boltzmann populaton n_(r,s)/N for any level r,s of E,
% input: temperature, boltzmann's constant, electron pressure, ion stage, ionization level

saha = saha_E(temp,p_e);
sahabolt = saha(1,ion) * boltz_E(temp,k,ion,lvl); 
    %using a designated vector element of saha according to ion stage

end

