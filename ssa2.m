%% SSA chapter 2

clc
close all
format long                             % prints double precision

% Constants and variables
k = 8.61734e-5;                         % Boltzmann's constant [eV/K]
kerg = 1.380658e-16;                    % Boltzmann's constant [erg K]
h = 6.62607e-27;                        % Planck's constant [erg s]
m_e = 9.109390e-28;                     % Electron mass [g]

p_e = 1e3;                              % electron pressure [dyne/cm^2]
temp = [5000 10000 20000];              % array of temperatures [K]


%% Boltzmann distribution

rstage = 1;                     % constant at ground state r = 1

for t = temp                    % running both functions for temperatures = 5000, 10 000, 20 000
u = partfunc_E(t,k);            % calling partitioning function with the input arguments t,k

disp(['T = ', num2str(t), ' K:', newline])
disp('U_r')     % header
disp(u.')       % display partition functions
disp('n_rs')    % display boltzmann

for s = 1:10                        % running for levels s = 1,..,10
    relnrs = boltz_E(t,k,rstage,s); % run Boltzmann function for levels s = 1...10
    disp(['    ', num2str(relnrs)])
end

disp(newline)
end

%% Saha distribution

for t = temp
disp(['T = ', num2str(t), ' K:', newline])
n_frac = saha_E(t,p_e);
disp(n_frac.')
end


%% Saha-Boltzmann

sb = zeros([1 5]);      % create a 5 zero-element array

for t = temp           % run for temperatures = 5000, 10 000, 20 000
 disp(['T = ', num2str(t), ' K:', newline])
 for s = 1:5        % run from ionization level 1-5
    sb = sahabolt_E(t,k,p_e,1,s);       % ion stage = ground stage = 1
    disp(sb)        % checking numerical values to see that code works
 end
end


%% Payne curves

ion = 1; 
p_ep = 131;                % Payne's electron pressure [dyne/cm^2]
temps = (0:1000:30000);    % length = 31
pop = zeros([5 31]);       % 5x31 matrix of zero-elements


 for T = 1:30           % goes through the array T
     for r = 1:4
        pop(r,T) = sahabolt_E(temps(T),k,p_ep,r,1);
     end
 end


% plot

figure
semilogy(temps,pop(1,:),'Color',[50/255 205/255 50/255],'Displayname', 'E^{+}')

hold on
grid on
grid minor

semilogy(temps,pop(2,:),'Color',[255/255 140/255 0/255],'Displayname', 'E^{2+}')
semilogy(temps,pop(3,:),'Color',[70/255 130/255 180/255],'Displayname', 'E^{3+}')
semilogy(temps,pop(4,:),'Color',[199/255 21/255 133/255],'Displayname', 'E^{4+}')

xlabel('Temperature [K]','Fontsize',12)
ylabel('Population','Fontsize',12)     % Stellar line strength
ylim([10e-3 1])

legend show
lgd.FontSize = 14;


%% Adding higher ionization levels

pop2 = zeros([5 31]);
pop3 = zeros([5 31]);
pop4 = zeros([5 31]);

 for T = 1:30
     for r = 1:4
        pop2(r,T) = sahabolt_E(temps(T),k,p_ep,r,2);
     end
 end

 for T = 1:30
     for r = 1:4
        pop3(r,T) = sahabolt_E(temps(T),k,p_ep,r,3);
     end
 end

 for T = 1:30
     for r = 1:4
        pop4(r,T) = sahabolt_E(temps(T),k,p_ep,r,4);
     end
 end


% plot

figure

p1 = semilogy(temps,pop(1,:),'Color',[50/255 205/255 50/255],'Displayname', 's = 1');

hold on
grid on
grid minor

semilogy(temps,pop(2:4,:),'Color',[50/255 205/255 50/255])

p2 = semilogy(temps,pop2(1,:),'Color',[255/255 140/255 0/255],'Displayname', 's = 2');
semilogy(temps,pop2(2:4,:),'Color',[255/255 140/255 0/255])

p3 = semilogy(temps,pop3(1,:),'Color',[70/255 130/255 180/255],'Displayname', 's = 3');
semilogy(temps,pop3(2:4,:),'Color',[70/255 130/255 180/255])

p4 = semilogy(temps,pop4(1,:),'Color',[199/255 21/255 133/255],'Displayname', 's = 4');
semilogy(temps,pop4(2:4,:),'Color',[199/255 21/255 133/255])


xlabel('Temperature [K]','Fontsize',12)
ylabel('Population','Fontsize',12)
ylim([10e-3 1])

leg = legend([p1 p2 p3 p4]);
lgd.FontSize = 14;
title(leg,'Ionization level')

 
 
%% Hydrogen Saha-Boltzmann

sbH = sahabolt_H(5000,1e2,1);
 

%% Solar Ca^+ K vs H alpha line strength

ionCaH = 1; 
p_CaH = 1e2;                    % Electron pressure [dyne/cm^2]
rstage = 1;                     % constant at ground state r = 1
tempCaH = (1000:100:20000);
sbCaH = zeros([1 5]);                       % 5 zero-element array

popCaH = zeros([5 length(tempCaH)]);        % 5x191 matrix of zero-elements
popCaH2 = zeros([5 length(tempCaH)]);
popCaH3 = zeros([5 length(tempCaH)]);
popCaH4 = zeros([5 length(tempCaH)]);

% Saha-Boltzmann
for t = tempCaH                         % running both functions for temperatures = 5000, 10 000, 20 000
    u = partfunc_CaH(t,k);              % calling partitioning function with the input arguments t,k

    for s = 1:10                            % running for levels s = 1,..,10
        relnrs = boltz_CaH(t,k,rstage,s);   % run Boltzmann function for levels s = 1...10
    end

n_frac = saha_CaH(t,p_e);

    for s2 = 1:5        % run from ionization level 1-5
        sbCaH = sahabolt_CaH(t,k,p_e,1,s2);       % ion stage = ground stage = 1
    end
end

 for T2 = 1:190         % goes through the array T-1
     for r2 = 1:4
        popCaH(r2,T2) = sahabolt_CaH(tempCaH(T2),k,p_CaH,r2,1);
        popCaH2(r2,T2) = sahabolt_CaH(tempCaH(T2),k,p_CaH,r2,2);
        popCaH3(r2,T2) = sahabolt_CaH(tempCaH(T2),k,p_CaH,r2,3);
        popCaH4(r2,T2) = sahabolt_CaH(tempCaH(T2),k,p_CaH,r2,4);
     end
 end

% Find strength ratio
CaH = zeros([1 length(tempCaH)]);
Caabund = 2e-6;

for i = 1:190
    NCa = sahabolt_CaH(tempCaH(i),k,p_CaH,2,1);   % is equal to sahabolt_Ca
    NH = sahabolt_H(tempCaH(i),p_CaH,2);
    CaH(i) = NCa.*Caabund./NH;
end

% Estimate solar line strength
loc = find(tempCaH == 5000);
ratio = CaH(loc);
disp(['Ca/H ratio at 5000 K = ',num2str(ratio)])



%% Ca plots

figure
semilogy(tempCaH,popCaH(1,:),'Color',[50/255 205/255 50/255],'Displayname', 'Ca^{+}')

hold on
grid on
grid minor

semilogy(tempCaH,popCaH(2,:),'Color',[255/255 140/255 0/255],'Displayname', 'Ca^{2+}')
semilogy(tempCaH,popCaH(3,:),'Color',[70/255 130/255 180/255],'Displayname', 'Ca^{3+}')
semilogy(tempCaH,popCaH(4,:),'Color',[199/255 21/255 133/255],'Displayname', 'Ca^{4+}')

xlabel('Temperature [K]','Fontsize',12)
ylabel('Population','Fontsize',12)     % Stellar line strength
ylim([10e-11 1])
legend show
lgd.FontSize = 14;

% all ionization energies
figure

cah1 = semilogy(tempCaH,popCaH(1,:),'Color',[50/255 205/255 50/255],'Displayname', 's = 1');

hold on
grid on
grid minor

   semilogy(tempCaH,popCaH(2:4,:),'Color',[50/255 205/255 50/255])

cah2 = semilogy(tempCaH,popCaH2(1,:),'Color',[255/255 140/255 0/255],'Displayname', 's = 2');
   semilogy(tempCaH,popCaH2(2:4,:),'Color',[255/255 140/255 0/255])

cah3 = semilogy(tempCaH,popCaH3(1,:),'Color',[70/255 130/255 180/255],'Displayname', 's = 3');
   semilogy(tempCaH,popCaH3(2:4,:),'Color',[70/255 130/255 180/255])

cah4 = semilogy(tempCaH,popCaH4(1,:),'Color',[199/255 21/255 133/255],'Displayname', 's = 4');
   semilogy(tempCaH,popCaH4(2:4,:),'Color',[199/255 21/255 133/255])


xlabel('Temperature [K]','Fontsize',12)
ylabel('Population','Fontsize',12)
ylim([10e-11 1])

leg = legend([cah1 cah2 cah3 cah4]);
lgd.FontSize = 14;
title(leg,'Ionization level')

%% Plot CaII K/ H alpha strength ratio
figure
semilogy(tempCaH,CaH,'Color',[50/255 205/255 50/255])
hold on
grid on
grid minor

xlabel('Temperature [K]','Fontsize',12)
ylabel('Ca II K / H \alpha','Fontsize',12)

legend('Strength ratio Ca II K / H \alpha')

%% Temperature sensitivity

tempc = (2000:100:12000);
dNCadT = zeros([1 length(tempc)]);
dNHdT = zeros([1 length(tempc)]);
dT = 1;

for i = 1:100
    NCa = sahabolt_CaH(tempc(i),k,p_CaH,2,1);
    NCa2 = sahabolt_CaH(tempc(i)-dT,k,p_CaH,2,1);
    
    dNCadT(i) = (NCa - NCa2)./(dT.*NCa);
    NH = sahabolt_H(tempc(i),1e2,2);
    NH2 = sahabolt_H(tempc(i)-dT,1e2,2);
    dNHdT(i) = (NH-NH2)./(dT.*NH);
end

% plot

figure
semilogy(tempc,abs(dNHdT),'Color',[50/255 205/255 50/255],'Displayname','H')

hold on
grid on
grid minor

semilogy(tempc,abs(dNCadT),'Color',[255/255 140/255 0/255],'Displayname','Ca^+ K')

xlabel('Temperature [K]','Fontsize',12)
ylabel('| (\Delta n_{(r,s)} / \Delta T) |','Fontsize',12)

legend show

%%

NCa = zeros([1 length(tempc)]);
NH = zeros([1 length(tempc)]);

for i = 1:100
    NCa(i) = sahabolt_CaH(tempc(i),k,p_CaH,2,1);    % Ca ion ground state
    NH(i) = sahabolt_H(tempc(i),p_CaH,2);           % H atom 2nd level
end

figure
hold on
grid on
grid minor

plot(tempc,NH./max(NH),'Color',[50/255 205/255 50/255],'Displayname','H');
plot(tempc,NCa./max(NCa),'Color',[255/255 140/255 0/255],'Displayname','Ca^+');

xlabel('Temperature [K]','Fontsize',12)
ylabel('Relative population','Fontsize',12)

legend show

%% Hot vs cool stars

for T3 = 2000:2000:20000
    sbHC = sahabolt_H(T3,1e2,1);
    disp([num2str(T3),'   ',num2str(sbHC)])
end

%% with plot

temphc = (1e3:1e2:2e4);
nH = zeros([1 length(temphc)]);

for i = 1:190
    nH(i) = sahabolt_H(temphc(i),1e2,1);
end

figure
hold on
grid on
grid minor

plot(temphc,nH,'Color',[50/255 205/255 50/255],'Displayname','H')
plot(temphc(83)-0.006422110050317,nH(83)-0.006422110050317,'x','Displayname','T = 9200 K')
disp(['T = ',num2str(temphc(83)-0.006422110050317),' K'])

xlabel('Temperature [K]','Fontsize',12)
ylabel('Neutral hydrogen fraction','Fontsize',12)

legend show



