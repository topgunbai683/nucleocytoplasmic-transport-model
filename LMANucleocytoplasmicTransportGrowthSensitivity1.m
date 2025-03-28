% function LMANucleocytoplasmicTransportGrowthSensitivity()
% Simulate a simplified model of Ran-mediated nucleocytoplasmic transport
% coupled with translation and growth
% using the law of mass action approximation.
% Based on: Wang et al. Thermodynamic Paradigm for Solution Demixing 
% Inspired by Nuclear Transport in Living Cells
% and our previous work.

% Sensitivity analysis using Pearson, Spearman, and PRCC.
% pv - parameters that are variable and whose sensitivities are analyzed.
% pf, phif - parameters that are assumed to be fixed and whose sensitivities 
% are not analyzed at the moment.
% Based on: N. G. Cogan. Mathematical Modeling the Life Sciences.

% X. B. 2024-10-05 Initial work

%% Steady-state solution of "pseudo-concentrations" in nmol*L^(-1)
% The following steady-state solution gives a roughly 8% NCcell.
% Where p.k2f0 = 5e-3*pf(4)^(-1) and pv(7)0 = 5e-2*pf(4), pf(4) = 1e2.
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
zsn = zsteady(4)+zsteady(5)+zsteady(11)+zsteady(12);
zst = zsteady(4)+zsteady(6)+zsteady(7)+zsteady(11)+zsteady(13)+zsteady(14);
zsran = zsteady(2)+zsteady(3)+zsteady(6)+zsteady(9)+zsteady(10)+zsteady(13);
zsntf = zsteady(1)+zsteady(3)+zsteady(8)+zsteady(10);
zstotal = zcy+zsn+zst+zsran+zsntf;

%% Parameters
% Vectors to store nominal values
pv = zeros(1,17); % Variable
pf = zeros(1,4); % Fixed
phif = zeros(1,2); % Fixed

% Gene fractions
% Ribosomes; taken from previous growth model
phif(1) = 0.2; % phir
% Nuclear proteins (that is, cargo proteins)
phif(2) = zsn/zstotal; % phin
% phin = 0.0735;
% phin = 0.074;
% phin = 0.075;
% phin = 0.08;
% phin = 0.085;
% phin = 0.2;
% Importins
pv(15) = zst/zstotal; % phit
% Ran
pv(16) = zsran/zstotal; % phiran
% NTF2
pv(17) = zsntf/zstotal; % phintf
% Other cytoplasmic proteins
phic = 1-(phif(1)+phif(2)+pv(15)+pv(16)+pv(17));

% Time
ts = 0;
% te = 1e1;
% te = 5e1;
% te = 1e2;
% te = 5e2;
% te = 1e3;
% te = 2e3;
te = 5e3;
% te = 6e3;
% te = 1e4;
% te = 1e5;
% te = 1e6;
tspan = [ts, te]; % Time interval in seconds
% Avogadro constant
N_A = 6.02214076e23; % mol^(-1)
% Reference volume Vref and reference number of molecules Nref
% Vref is also the initial cytoplasmic volume.
pf(3) = 500*(1e-6)^3*1e3; % L
% Determination of the "dry volume" V0 in the cytoplasmic volume
% Note here that V0 should not dependent on zcy.
% pf(1) = 1e-1; % "Initial dry volume fraction", tentative
pf(1) = 1e-2; % "Initial dry volume fraction", tentative
V0 = pf(1)*pf(3);
% Initial NPC permeability
% a0 = 5e1*(1e-6)^3*1e3; % L*s^(-1)
a0 = 1e2*(1e-6)^3*1e3; % L*s^(-1)
% a0 = 1.5e2*(1e-6)^3*1e3; % L*s^(-1)
% a0 = 1e3*(1e-6)^3*1e3; % L*s^(-1)
% GTP:GDP ratio theta
% pf(4) = 5e1;
pf(4) = 1e2;
% pf(4) = 1e3;
% Translation rate
% This may be calculated from the doubling time.
% pv(13) = 5e-4;
% pv(13) = 1e-3;
% pv(13) = 4e-3;
pv(13) = 5e-3;
% pv(13) = 1e-2;
% Reaction rates
pv(1) = 0.1; % nmol^(-1)*L*s^(-1)
pv(2) = 2.5; % s^(-1)
% pv(3) = 5e1; % s^(-1)
% pv(3) = 1e2; % s^(-1)
pv(3) = 3e2; % s^(-1)
pv(4) = 1; % mol^(-1)*L*s^(-1)
pv(5) = 0.1; % nmol^(-1)*L*s^(-1)
% pv(6) = 1; % s^(-1)
pv(6) = 1e-1; % s^(-1)
% pv(6) = 1e-2; % s^(-1)
% pv(6) = 1e-6; % s^(-1)
% pv(7) = 50; % s^(-1)
% pv(7) = 1e2; % s^(-1)
pv(7) = 2e2; % s^(-1)
pv(8) = 1; % nmol^(-1)*L*s^(-1)
pv(9) = 20; % s^(-1)
pv(10) = 1; % nmol^(-1)*L*s^(-1)
pv(11) = 1; % nmol^(-1)*L*s^(-1)
pv(12) = 20; % s^(-1)

%%  Initial states of the "pseudo-concentrations"
z0 = [zsteady; phif(1)*zstotal; phic*zstotal];
% Test whether the NC ratio homeostasis is maintained.
% z0p means "z0 perturbed".
z0p = z0; % No perturbation in the initial NC ratio
% z0p = [zsteady; phir*zstotal; 1.1*phic*zstotal];
% z0p = [zsteady; phir*zstotal; 1.2*phic*zstotal];
% z0p = [zsteady; phir*zstotal; 0.9*phic*zstotal];
% z0p = [zsteady; phir*zstotal; 0.8*phic*zstotal];
% Test the effect of the order of magnitude
% z0 = 1e-1*z0;
% Determination of the parameter Ccy in the cytoplasmic volume
zcyt0 = sum(z0(1:7))+sum(z0(15:16));
pf(2) = (pf(3)-V0)/zcyt0;
% Determination of the parameter Cnp in the NPC permeability
NCcyto0 = sum(z0(8:14))/zcyt0;
Vn0 = NCcyto0*pf(3);
% Scaling with respect to nuclear surface area
% pv(14) = a0*Vn0^(-2/3);
% Scaling with respect to nuclear volume
pv(14) = a0*Vn0^(-1);
% Scaling with respect to cytoplasmic protein number
% pv(14) = a0*zcyt0^(-1);

%% Setup for sensitivity analysis
set(groot,'DefaultAxesTickLabelInterpreter','latex');
set(groot,'DefaultLegendInterpreter','latex');
%{
Params_names={'$k_1^+$','$k_1^-$','$k_{2,0}^+$','$k_2^-$','$k_4^+$',...
              '$k_4^-$','$k_{5,0}^+$','$k_5^-$','$k_7^+$','$k_7^-$',...
              '$k_{10}^+$', '$k_{10}^-$','$k_t$','$C_{np}$','$\phi_{Imp}$',...
              '$\phi_{Ran}$','$\phi_{NTF}$'};
%}

%% Generate parameters using Latin hypercube sample
% Num_samples = 10; % For testing the codes
% Num_samples = 50;
% Num_samples = 500;
Num_samples = 1500;

% %{
% Uniform distribution
punif = lhsdesign(Num_samples,length(pv));
% l_bounds = .95*pv;
% u_bounds = 1.05*pv;
l_bounds = .9*pv;
u_bounds = 1.1*pv;
pvm = zeros(size(punif));
for k = 1:Num_samples
    pvm(k,:) = (u_bounds-l_bounds).*punif(k,:)+l_bounds;
end
% %}

% Normal distribution
% pvm = abs(lhsnorm(pv,diag(.1*pv),Num_samples));

%% Evaluate the model for each parameter set
NCcellv = zeros(Num_samples,1); % Column vector
for k = 1:Num_samples
    pvk = pvm(k,:);
    [t,z] = ode15s(@(t,z) odefun_nucleocytoplasmictransportgrowthsensitivity1(t,z,pvk,pf,phif), tspan, z0p);
    NCcell = sum(z(:,8:14),2)./sum(z,2); % This could be simplified
    NCcellv(k) = NCcell(end);
end

%% Sensitivity analysis
ccp = zeros(1,length(pv));
ccs = zeros(1,length(pv));
for k = 1:length(pv)
    CCp=corr(pvm(:,k),NCcellv,'Type','Pearson');
    CCs=corr(pvm(:,k),NCcellv,'Type','Spearman');
    ccp(k)=CCp(1,end);
    ccs(k)=CCs(1,end);
end
correlation_matrix=[pvm,NCcellv];
PRCC=partialcorr(correlation_matrix);
prcc = PRCC(:,end);

%% Plots
figure(1)
bar_width=.35;
bar([1:length(pv)]-0*bar_width,ccp,bar_width,'c')
xticks(1:length(pv))
xticklabels({'$k_1^+$','$k_1^-$','$k_{2,0}^+$','$k_2^-$','$k_4^+$',...
              '$k_4^-$','$k_{5,0}^+$','$k_5^-$','$k_7^+$','$k_7^-$',...
              '$k_{10}^+$', '$k_{10}^-$','$k_t$','$C_{np}$','$\phi_{Imp}$',...
              '$\phi_{Ran}$','$\phi_{NTF}$'})
title('Pearson')

figure(2)
bar_width=.35;
bar([1:length(pv)]-0*bar_width,ccs,bar_width,'c')
xticks(1:length(pv))
xticklabels({'$k_1^+$','$k_1^-$','$k_{2,0}^+$','$k_2^-$','$k_4^+$',...
              '$k_4^-$','$k_{5,0}^+$','$k_5^-$','$k_7^+$','$k_7^-$',...
              '$k_{10}^+$', '$k_{10}^-$','$k_t$','$C_{np}$','$\phi_{Imp}$',...
              '$\phi_{Ran}$','$\phi_{NTF}$'})
title('Spearman')

figure(3)
bar_width=.35;
bar([1:length(pv)]-0*bar_width,prcc(1:length(pv)),bar_width,'c')
xticks(1:length(pv))
xticklabels({'$k_1^+$','$k_1^-$','$k_{2,0}^+$','$k_2^-$','$k_4^+$',...
              '$k_4^-$','$k_{5,0}^+$','$k_5^-$','$k_7^+$','$k_7^-$',...
              '$k_{10}^+$', '$k_{10}^-$','$k_t$','$C_{np}$','$\phi_{Imp}$',...
              '$\phi_{Ran}$','$\phi_{NTF}$'})
title('PRCC')
