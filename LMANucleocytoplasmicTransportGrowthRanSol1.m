% function LMANucleocytoplasmicTransportGrowthRanSol1()
% Simulate a simplified model of Ran-mediated nucleocytoplasmic transport
% coupled with translation and growth
% using the law of mass action approximation.
% Based on: Wang et al. Thermodynamic Paradigm for Solution Demixing 
% Inspired by Nuclear Transport in Living Cells
% and our previous work.

% Cytoplasmic volume Vcy is set to be an affine function of the total number 
% of cytoplasmic proteins.
% Nuclear volume Vn is variable and is determined by the osmotic balance
% at the nuclear envelope.
% The variables z(i) are "pseudo-concentrations" measured in nM = nmol*L^(-1).
% That is, z(i)=(x(i)/(N_A*Vref))*1e9, where x(i) are the numbers of molecules.
% For cytoplasmic species the true concentration is z(i)*(Vref/Vcy).
% For nuclear species the true concentration is z(i)*(Vref/Vn).

% For now gene fractions is determined from steady-state solutions obtained
% from LMANucleocytoplasmicTransportNCVSol1.
% Therefore these gene fractions may not be biologically realistic.
% Later we may use yeast genome data to determine the gene fractions.

% Assume Cytoplasmic RanGDP can diffuse through NPCs with a much lower
% diffusion rate than others.

% X. B. 2025-02-20 Initial work

%% Steady-state solution of "pseudo-concentrations" in nmol*L^(-1)
% Make sure the reference volume Vref of the steady-state solution 
% is the same as the p.Vref in the Parameters section.
%{
% The following steady-state solution gives a roughly 8% NCcell.
% Where p.k2f = 50*p.theta^(-1) and p.k5f = 50*p.theta, p.theta = 100.
% Some values could be modified away from the original steady state.
zcy = 4.5e+04; % All the other cytoplasmic proteins
zsteady = [5.972011504852738e+02; % z(1):Cytoplasmic NTF2
           13.531213929293655;    % z(2):Cytoplasmic RanGDP
           3.228230883833524e+02; % z(3):Cytoplasmic NTF2-RanGDP complex
           2.189260911333862e+02; % z(4):Cytoplasmic Importin-cargo complex
           8.357013161413893e+02; % z(5):Cytoplasmic cargo
           1.438442221275832e+02; % z(6):Cytoplasmic Importin-RanGTP complex
           5.239337031795179;     % z(7):Cytoplasmic Importin
           % 5e1*5.239337031795179;     % z(7):Cytoplasmic Importin
           % 1e2*5.239337031795179;     % z(7):Cytoplasmic Importin
           52.360100015880676;    % z(8):Nuclear NTF2
           2.292349818909228e+02; % z(9):Nuclear RanGTP
           27.615658776785460;    % z(10):Nuclear NTF2-RanGDP complex
           19.030750403370470;    % z(11):Nuclear Importin-cargo complex
           3.756341809916630e+03; % z(12):Nuclear cargo
           % 1.1*3.756341809916630e+03; % z(12):Nuclear cargo
           % 0.9*3.756341809916630e+03; % z(12):Nuclear cargo
           12.950741127214199;    % z(13):Nuclear Importin-RanGTP complex
           0.008808052096420      % z(14):Nuclear Importin
          ];
%}
%{
% The following steady-state solution gives a roughly 8% NCcell.
% Where p.k2f = 3e2*p.theta^(-1) and p.k5f = 2e2*p.theta, p.theta = 100.
zcy = 4.5e+04; % All the other cytoplasmic proteins
zsteady = [7.084633293094366e+02; % z(1):Cytoplasmic NTF2
           7.538876888059733;     % z(2):Cytoplasmic RanGDP
           2.115627868991066e+02; % z(3):Cytoplasmic NTF2-RanGDP complex
           2.745527128293247e+02; % z(4):Cytoplasmic Importin-cargo complex
           2.102982899277234e+02; % z(5):Cytoplasmic cargo
           67.346946203161070;    % z(6):Cytoplasmic Importin-RanGTP complex
           26.110789036010832;    % z(7):Cytoplasmic Importin
           % 5*26.110789036010832;    % z(7):Cytoplasmic Importin
           63.841448905170540;    % z(8):Nuclear NTF2
           4.393065992379086e+02; % z(9):Nuclear RanGTP
           16.132419627200512;    % z(10):Nuclear NTF2-RanGDP complex
           23.865673120591143;    % z(11):Nuclear Importin-cargo complex
           3.491283305980321e+03; % z(12):Nuclear cargo
           8.111990323557386;     % z(13):Nuclear Importin-RanGTP complex
           0.011884106652177      % z(14):Nuclear Importin
           ];
%}
% 
%{
% The following steady-state solution gives a roughly 8% NCcell.
% Where p.k2f = 3e2*p.theta^(-1) and p.k5f = 2e2*p.theta, p.theta = 100.
zcy = 4.5e5; % All the other cytoplasmic proteins
zsteady = [2.846714874736455e+02; % z(1):Cytoplasmic NTF2
           56.319611851871440;    % z(2):Cytoplasmic RanGDP
           6.352681463875213e+02; % z(3):Cytoplasmic NTF2-RanGDP complex
           2.180934549708345e+03; % z(4):Cytoplasmic Importin-cargo complex
           5.773781741439327e+02; % z(5):Cytoplasmic cargo
           1.423272202896009e+03; % z(6):Cytoplasmic Importin-RanGTP complex
           75.546116323837310;    % z(7):Cytoplasmic Importin
           31.339914321180680;    % z(8):Nuclear NTF2
           2.705795700162941e+03; % z(9):Nuclear RanGTP
           48.720451633213834;    % z(10):Nuclear NTF2-RanGDP complex
           1.898020390965505e+02; % z(11):Nuclear Importin-cargo complex
           3.651188492576560e+04; % z(12):Nuclear cargo
           1.304298835280764e+02; % z(13):Nuclear Importin-RanGTP complex
           0.009048047175478      % z(14):Nuclear Importin
           ];
%}
% %{
% The following steady-state solution gives a roughly 8% NCcell.
% Where p.k2f0 = 5e-3*p.theta^(-1) and p.k5f0 = 5e-2*p.theta, p.theta = 1e2.
zcy = 4.5e+04; % All the other cytoplasmic proteins
zsteady = [7.064840782810440e+02; % z(1):Cytoplasmic NTF2
           7.603564787986136;     % z(2):Cytoplasmic RanGDP
           2.132206262357328e+02; % z(3):Cytoplasmic NTF2-RanGDP complex
           2.537841863194739e+02; % z(4):Cytoplasmic Importin-cargo complex
           2.444211453520565e+02; % z(5):Cytoplasmic cargo
           93.331557026060680;    % z(6):Cytoplasmic Importin-RanGTP complex
           20.766140018615374;    % z(7):Cytoplasmic Importin
           63.482010596670170;    % z(8):Nuclear NTF2
           4.090801736656363e+02; % z(9):Nuclear RanGTP
           16.813269226681770;    % z(10):Nuclear NTF2-RanGDP complex
           22.156755477852037;    % z(11):Nuclear Importin-cargo complex
           % 3.538354588658927e+03; % z(12):Nuclear cargo
           3.6e+03; % z(12):Nuclear cargo
           % 4e+03; % z(12):Nuclear cargo
           % 3e+03; % z(12):Nuclear cargo
           9.950422654274737;     % z(13):Nuclear Importin-RanGTP complex
           0.010933933187063      % z(14):Nuclear Importin
           ];
% %}
zsn = zsteady(4)+zsteady(5)+zsteady(11)+zsteady(12);
zst = zsteady(4)+zsteady(6)+zsteady(7)+zsteady(11)+zsteady(13)+zsteady(14);
zsran = zsteady(2)+zsteady(3)+zsteady(6)+zsteady(9)+zsteady(10)+zsteady(13);
zsntf = zsteady(1)+zsteady(3)+zsteady(8)+zsteady(10);
zstotal = zcy+zsn+zst+zsran+zsntf;
% Modify zsteady after the gene fractions are calculated
% zsteady(12) = 3.7e+04; % z(12):Nuclear cargo
% Gene fractions
% Ribosomes; taken from previous growth model
phir = 0.2; 
% Nuclear proteins (that is, cargo proteins)
phin = zsn/zstotal;
% phin = 0.0735;
% phin = 0.074;
% phin = 0.075;
% phin = 0.08;
% phin = 0.085;
% phin = 0.2;
% Importins
phit = zst/zstotal;
% Ran
phiran = zsran/zstotal;
% phiran = 1.1*zsran/zstotal;
% phiran = 2*zsran/zstotal;
% phiran = 3*zsran/zstotal;
% NTF2
phintf = zsntf/zstotal;
% phintf = 1.1*zsntf/zstotal;
% phintf = 1e-4;
% phintf = 0; % Complete repression of NTF2 expression
% Other cytoplasmic proteins
phic = 1-(phir+phin+phit+phiran+phintf);

%% Parameters
% Time
ts = 0;
% te = 1e1;
% te = 5e1;
% te = 1e2;
% te = 5e2;
% te = 1e3;
% te = 2e3;
% te = 5e3;
% te = 6e3;
% te = 1e4;
te = 1.5e4;
% te = 2e4;
% te = 1e5;
% te = 1e5;
% te = 1e6;
tspan = [ts, te]; % Time interval in seconds
% Avogadro constant
p.N_A = 6.02214076e23; % mol^(-1)
% Reference volume Vref and reference number of molecules Nref
% Vref is also the initial cytoplasmic volume.
p.Vref = 500*(1e-6)^3*1e3; % L
p.Nref = p.N_A*p.Vref*1e-9; % nmol^(-1)*L
% Determination of the "dry volume" V0 in the cytoplasmic volume
% Note here that V0 should not dependent on zcy.
% nu = 1e-1; % "Initial dry volume fraction", tentative
% p.V0 = nu*p.Vref;
% p.nu = 1e-1; % "Initial dry volume fraction", tentative
p.nu = 1e-2; % "Initial dry volume fraction", tentative
p.V0 = p.nu*p.Vref;
% Initial NPC permeability
% p.a0 = 5e1*(1e-6)^3*1e3; % L*s^(-1)
p.a0 = 1e2*(1e-6)^3*1e3; % L*s^(-1)
% p.a0 = 1.5e2*(1e-6)^3*1e3; % L*s^(-1)
% p.a0 = 1e3*(1e-6)^3*1e3; % L*s^(-1)
% Scaling factor for NPC permeability of RanGDP
% p.Cran = 0;
% p.Cran = 1e-2;
p.Cran = 3e-2;
% GTP:GDP ratio
% p.theta = 5e1;
p.theta = 1e2;
% p.theta = 1e3;
% Translation rate
% This may be calculated from the doubling time.
% p.kt = 5e-4;
% p.kt = 1e-3;
% p.kt = 4e-3;
p.kt = 5e-3;
% p.kt = 1.1*5e-3;
% p.kt = 1e-2;
% Reaction rates
p.k1f = 0.1; % nmol^(-1)*L*s^(-1)
p.k1b = 2.5; % s^(-1)
% p.k2f = 50*p.theta^(-1); % s^(-1)
% p.k2f = 1e2*p.theta^(-1); % s^(-1)
p.k2f = 3e2*p.theta^(-1); % s^(-1)
% p.k2f = 4.5e2*p.theta^(-1); % s^(-1)
% p.k2f = 6e2*p.theta^(-1); % s^(-1)
% p.k2f = 4.5e3*p.theta^(-1); % s^(-1)
p.k2b = 1; % mol^(-1)*L*s^(-1)
p.k4f = 0.1; % nmol^(-1)*L*s^(-1)
% p.k4b = 1; % s^(-1)
p.k4b = 1e-1; % s^(-1)
% p.k4b = 1e-2; % s^(-1)
% p.k4b = 1e-6; % s^(-1)
% p.k5f = 50*p.theta; % s^(-1)
% p.k5f = 1e2*p.theta; % s^(-1)
p.k5f = 2e2*p.theta; % s^(-1)
% p.k5f = 3e2*p.theta; % s^(-1)
% p.k5f = 2e3*p.theta; % s^(-1)
p.k5b = 1; % nmol^(-1)*L*s^(-1)
p.k7f = 20; % s^(-1)
p.k7b = 1; % nmol^(-1)*L*s^(-1)
p.k10f = 1; % nmol^(-1)*L*s^(-1)
p.k10b = 20; % s^(-1)
p.k12f = p.kt*phir; % s^(-1)
p.k13f = p.kt*phin; % s^(-1)
p.k14f = p.kt*phic; % s^(-1)
p.k15f = p.kt*phit; % s^(-1)
p.k16f = p.kt*phiran; % s^(-1)
p.k17f = p.kt*phintf; % s^(-1)
p.k19f = 1e2*p.theta; % s^(-1)
p.k19b = 1; % s^(-1)

%%  Initial states of the "pseudo-concentrations"
z0 = [zsteady; phir*zstotal; phic*zstotal; 0];
% Test whether the NC ratio homeostasis is maintained.
% z0p means "z0 perturbed".
z0p = z0; % No perturbation in the initial NC ratio
% z0p = [zsteady; phir*zstotal; 1.1*phic*zstotal; 0];
% z0p = [zsteady; phir*zstotal; 1.2*phic*zstotal; 0];
% z0p = [zsteady; phir*zstotal; 0.9*phic*zstotal; 0];
% z0p = [zsteady; phir*zstotal; 0.8*phic*zstotal; 0];
% Test the effect of the order of magnitude
% z0 = 1e-1*z0;
% Determination of the parameter Ccy in the cytoplasmic volume
zcyt0 = sum(z0(1:7))+sum(z0(15:16));
p.Ccy = (p.Vref-p.V0)/zcyt0;
% Determination of the parameter Cnp in the NPC permeability
NCcyto0 = (sum(z0(8:14))+z0(17))/zcyt0;
Vn0 = NCcyto0*p.Vref;
% Scaling with respect to nuclear surface area
% p.Cnp = p.a0*Vn0^(-2/3);
% Scaling with respect to nuclear volume
p.Cnp = p.a0*Vn0^(-1);
% p.Cnp = 1.1*p.Cnp;
% Scaling with respect to cytoplasmic protein number
% p.Cnp = p.a0*zcyt0^(-1);
% Set initial NTF2 concentration to 0
% z0p(1) = 0; z0p(3) = 0; z0p(8) = 0; z0p(10) = 0;

% Experimenting with scaling
% p.Cnp = p.a0*Vn0^(-4/3);
% p.Cnp = p.a0*Vn0^(-0.9);
% p.Cnp = p.a0*Vn0^(-1.1);

%% Solve the ODE
% options=odeset('RelTol',1e-6,'AbsTol',1e-6,'Stats','on');
% [t,z] = ode45(@(t,z) odefun_nucleocytoplasmictransportgrowthran1(t,z,p), tspan, z0p);
[t,z] = ode15s(@(t,z) odefun_nucleocytoplasmictransportgrowthran1(t,z,p), tspan, z0p);
% [t,z] = ode23s(@(t,z) odefun_nucleocytoplasmictransportgrowthran1(t,z,p), tspan, z0p);
% Calculate the N/C ratios
NCcyto = (sum(z(:,8:14),2)+z(:,17))./(sum(z(:,1:7),2)+sum(z(:,15:16),2));
NCcell = (sum(z(:,8:14),2)+z(:,17))./sum(z,2);
% Calculate total number of proteins
zn = z(:,4)+z(:,5)+z(:,11)+z(:,12);
zt = z(:,4)+z(:,6)+z(:,7)+z(:,11)+z(:,13)+z(:,14);
zran = z(:,2)+z(:,3)+z(:,6)+z(:,9)+z(:,10)+z(:,13)+z(:,17);
zntf = z(:,1)+z(:,3)+z(:,8)+z(:,10);
zcyt = sum(z(:,1:7),2)+sum(z(:,15:16),2);
Vcy = p.V0+p.Ccy*zcyt;
Vn = Vcy.*NCcyto;

%% Plot the results
%{
figure(1)
% plot(t,z(:,1),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,1),'Linewidth',2) % Number
plot(t,p.Vref*z(:,1)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Cytoplasmic NTF2') % Number
ylabel('Cytoplasmic NTF2 (nM)') % True concentration
%}
%{
figure(2)
% plot(t,z(:,2),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,2),'Linewidth',2) % Number
plot(t,p.Vref*z(:,2)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Cytoplasmic RanGDP') % Number
ylabel('Cytoplasmic RanGDP (nM)') % True concentration
%}
%{
figure(3)
% plot(t,z(:,3),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,3),'Linewidth',2) % Number
plot(t,p.Vref*z(:,3)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Cytoplasmic NTF2-RanGDP complex') % Number
ylabel('Cytoplasmic NTF2-RanGDP complex (nM)') % True concentration
%}
%{
figure(4)
% plot(t,z(:,4),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,4),'Linewidth',2) % Number
plot(t,p.Vref*z(:,4)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Cytoplasmic importin-cargo complex') % Number 
ylabel('Cytoplasmic importin-cargo complex (nM)') % True concentration
%}
%{
figure(5)
% plot(t,z(:,5),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,5),'Linewidth',2) % Number
plot(t,p.Vref*z(:,5)./Vcy,'Linewidth',2) % True concentration
% axis([0 inf 0 5e3])
xlabel('Time (s)')
% ylabel('Cytoplasmic cargo') % Number
ylabel('Cytoplasmic cargo (nM)') % True concentration
%}
%{
figure(6)
% plot(t,z(:,6),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,6),'Linewidth',2) % Number
plot(t,p.Vref*z(:,6)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Cytoplasmic importin-RanGTP complex') % Number
ylabel('Cytoplasmic importin-RanGTP complex (nM)') % True concentration
%}
%{
figure(7)
% plot(t,z(:,7),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,7),'Linewidth',2) % Number
xplot(t,p.Vref*z(:,7)./Vcy,'Linewidth',2) % True concentration
label('Time (s)')
% ylabel('Cytoplasmic importin') % Number
ylabel('Cytoplasmic importin (nM)') % True concentration
%}
%{
figure(8)
% plot(t,z(:,8),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,8),'Linewidth',2) % Number
plot(t,p.Vref*z(:,8)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear NTF2') % Number
ylabel('Nuclear NTF2 (nM)') % True concentration
%}
%{
figure(9)
% plot(t,z(:,9),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,9),'Linewidth',2) % Number
plot(t,p.Vref*z(:,9)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear RanGTP') % Number
ylabel('Nuclear RanGTP (nM)') % True concentration
%}
%{
figure(10)
% plot(t,z(:,10),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,10),'Linewidth',2) % Number
plot(t,p.Vref*z(:,10)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear NTF2-RanGDP complex') % Number 
ylabel('Nuclear NTF2-RanGDP complex (nM)') % True concentration
%}
%{
figure(11)
% plot(t,z(:,11),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,11),'Linewidth',2) % Number
plot(t,p.Vref*z(:,11)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear importin-cargo complex') % Number
ylabel('Nuclear importin-cargo complex (nM)') % True concentration
%}
%{
figure(12)
% plot(t,z(:,12),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,12),'Linewidth',2) % Number
plot(t,p.Vref*z(:,12)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear cargo') % Number
ylabel('Nuclear cargo (nM)') % True concentration
%}
%{
figure(13)
% plot(t,z(:,13),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,13),'Linewidth',2) % Number
plot(t,p.Vref*z(:,13)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear importin-RanGTP complex') % Number
ylabel('Nuclear importin-RanGTP complex (nM)') % True concentration
%}
%{
figure(14)
% plot(t,z(:,14),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,14),'Linewidth',2) % Number
plot(t,p.Vref*z(:,14)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear importin') % Number
ylabel('Nuclear importin (nM)') % True concentration
%}
%{
figure(15)
% plot(t,z(:,15),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,15),'Linewidth',2) % Number
plot(t,p.Vref*z(:,15)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear importin') % Number
ylabel('Ribosomes (nM)') % True concentration
%}
%{
figure(16)
% plot(t,z(:,16),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,16),'Linewidth',2) % Number
plot(t,p.Vref*z(:,16)./Vcy,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Other cytoplasmic proteins') % Number
ylabel('Other cytoplasmic proteins (nM)') % True concentration
%}
%{
figure(17)
% plot(t,z(:,17),'Linewidth',2)
% plot(t,p.N_A*p.Vref*1e-9*z(:,17),'Linewidth',2) % Number
plot(t,p.Vref*z(:,17)./Vn,'Linewidth',2) % True concentration
xlabel('Time (s)')
% ylabel('Nuclear RanGDP') % Number
ylabel('Nuclear RanGDP (nM)') % True concentration
%}
% %{
figure(18)
plot(t,NCcell,'Linewidth',2)
% axis([0 inf 0 0.2])
% axis([0 inf 0.05 0.2])
xlabel('Time (s)')
ylabel('Nuclear-to-cell ratio')
% %}
%{
figure(19)
plot(t,Vcy,'Linewidth',2)
xlabel('Time (s)')
ylabel('Cytoplasmic volume')
%}
%{
figure(20)
plot(t,Vn,'Linewidth',2)
xlabel('Time (s)')
ylabel('Nuclear volume')
%}

%% Save data
%{
save('NucleocytoplasmicTransport_AllNominal.mat','phir','phin',...
'phit','phiran','phintf','phic','p','t','z','NCcyto','NCcell','zn','zt',...
'zran','zntf','zcyt','Vcy','Vn')
%}
% %{
save(['NucleocytoplasmicTransport_phiran=' num2str(phiran) '_phintf=' num2str(phintf) '_Cran=' num2str(p.Cran) '.mat'],...
'phir','phin','phit','phiran','phintf','phic','p','t','z','NCcyto','NCcell',...
'zn','zt','zran','zntf','zcyt','Vcy','Vn')
% %}

