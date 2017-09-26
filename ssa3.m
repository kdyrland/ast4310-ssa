%% SSA chapter 3

clc
clearvars
close all
format long

wav = (1000:200:20800);
b = zeros([1 length(wav)]);

figure
hold on
grid on
grid minor

for T = 8000:(-200):5000
    b = planck(T,wav*1e-8);   %[K,cm]

    plot(wav*1e-8,b)

end


legend('show','T = 8000 K','T = 7800 K','T = 7600 K','T = 7400 K','T = 7200 K','T = 7000 K','T = 6800 K','T = 6600 K','T = 6400 K','T = 6200 K','T = 6000 K','T = 5800 K','T = 5600 K','T = 5400 K','T = 5200 K','T = 5000 K')

xlabel('Wavelength [cm]')
ylabel('B_{\lambda}(T)')

%% Emergent intensity

B = 2;
tau = (0.01:0.01:10);
intensity = zeros([1 length(tau)]);

figure
hold on
grid on
grid minor

for I0 = 4:-1:0
    intensity = I0 + exp(-tau) + B.*(1 - exp(-tau));
    plot(tau,intensity)
end

xlabel('Optical depth \tau')
ylabel('Intensity I')
legend('show','I(0) = 4','I(0) = 3','I(0) = 2','I(0) = 1','I(0) = 0')

%% Voigt profile
% NB: does not work properly,  so from here on i used python

clc
u = (-10:0.1:10);             % dimensionless units for wavelength scale
a = [0.001 0.01 0.1 1];
vau = zeros([length(a) length(u)]);

figure
hold on
grid on
grid minor

for i = 1:4
    vau(i,:) = voigt(a(i),abs(u));
    plot(u,vau(i,:))
end

xlabel('u')
ylabel('Voigt function')
ylim([0 1])
xlim([-10 10])
legend('a = 0.001','a = 0.01','a = 0.1','a = 1')



%% Schuster-Schwarzschild

Ts = 5700;          % solar surface temperature
Tl = 4200;          % solar T-min temperature = 'reversing layer'
a1 = 0.1;           % damping parameter
wav = 5000e-8;      % wavelength in cm
tau0 = 1;           % reversing layer thickness at line center
tau = zeros(size(u));
intensity = zeros(size(u));

figure
hold on
for i = 1:200
    tau = tau0 * voigt(a1, u(i));
    intensity(i) = planck(Ts,wav) .* exp(-tau) + planck(Tl,wav).*(1-exp(-tau));
end


% logtau0 = (-2:0.5:2);
% 
% for itau = 1:10
%     for i = 1:200
%         tau = 10^(logtau0(itau)) * voigt(a1, u(i));
%         intensity(i) = planck(Ts,wav) .* exp(-tau) + planck(Tl,wav).*(1-exp(-tau));
%     end
     plot(u,intensity)
     xlim([-15 15])
% end
% 



