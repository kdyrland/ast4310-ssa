function [u] = partfunc_CaH(temp,k)   % function [output_args] = func_name(input_args)
% partition functions Schadee element E
% input: temperature and boltzmann; k so we only need to define it in main script

    chiCaH = [6.113 11.871 50.91 67.15];  % Ca ionization energies [eV]
    u = zeros([1 4]);       % declare 4 zero-element array

    for rx = 1:4
         for s = 0:chiCaH(rx)-1     % (chi-1) bc r starts at 1, while s starts at 0
            u(1,rx) = u(1,rx) + exp(-s/(k*temp));
         end
    end

end