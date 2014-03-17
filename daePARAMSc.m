function params = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,factor)

% Parameters for rat signaling model
% params = ratPARAMS()
%
% Copyright 2011, Cardiac Systems Biology Lab, University of Virginia
%   JS: Jeff Saucerman  <jsaucerman@virginia.edu>
%   JY: Jason Yang      <jhyang@virginia.edu>
%   RA: Robert Amanfu   <rka2p@virginia.edu>
%
% Robert Amanfu
% 11/08/11

%% Signaling Parameters
disp('desensitization modification');
% Drug Concentrations
Ltot            = 0;              % (uM)ligand concentration
Atot            = 0;              % (uM)agonist concentration
FSK             = 0;              % (uM) forskolin concentration
IBMX            = 0;              % (uM) IBMX concentration

% ETC receptor module/Gs module
b1ARtot = 0.0132;
Gstot =  3.83;
kf_bARK = 1.1e-6;
kr_bARK = 2.2e-6;%2.2e-6;
k2_bARK = 5.5e-06;
Km_bARK = 6.26e-04;
kf_PKA =  3.6e-6;
kr_PKA = 2.2e-6;
k_G_act = 16e-3;
k_G_hyd = 0.8e-3;
k_G_reassoc = 1.21;
kf_Ra = 1;   % (1/[uM ms]) forward rate for switching between Ri and Ra
kr_Ra  = KR;   
kf_LRa1  = 1; 
kr_LRa1  = (alpha_L*KL);
kf_LRi = 1;   
kr_LRi = KL;
kf_LRa2 = 1;  
kr_LRa2 = (alpha_L*KR);
kf_RaG = 1; 
kr_RaG = KG;
kf_LRaG1  = 1;
kr_LRaG1 = (alpha_L*gamma_L*KL);
kf_LRaG2 = 1;
kr_LRaG2 = (gamma_L*KG);
kf_ARa1  = 1; 
kr_ARa1  = (alpha_A*KA);
kf_ARi = 1;   
kr_ARi = KA;
kf_ARa2 = 1;  
kr_ARa2 = (alpha_A*KR);
kf_ARaG1  = 1;
kr_ARaG1 = (alpha_A*gamma_A*KA);
kf_ARaG2 = 1;
kr_ARaG2 = (gamma_A*KG);

% AC module
ACtot           = 49.7e-3;       % (uM) total adenylyl cyclase
ATP             = 5e3;            % (uM) total ATP
k_AC_basal      = 0.2e-3;         % (1/ms) basal cAMP generation rate by AC
Km_AC_basal     = 1.03e3;         % (uM) basal AC affinity for ATP

Kd_AC_Gsa       =  0.2250;%0.4;            % (uM) Kd for AC association with Gsa
kf_AC_Gsa       = 1;              % (1/[uM ms]) forward rate for AC association with Gsa
kr_AC_Gsa       = Kd_AC_Gsa;      % (1/ms) reverse rate for AC association with Gsa

k_AC_Gsa        = 8.5e-3;         % (1/ms) basal cAMP generation rate by AC:Gsa
Km_AC_Gsa       = 315.0;          % (uM) AC:Gsa affinity for ATP

Kd_AC_FSK       = 44.0;           % (uM) Kd for FSK binding to AC
k_AC_FSK        = 7.3e-3;         % (1/ms) basal cAMP generation rate by AC:FSK
Km_AC_FSK       = 860.0;          % (uM) AC:FSK affinity for ATP

PDEtot          = 38.9e-3;       % (uM) total phosphodiesterase
k_cAMP_PDE      = 5e-3;           % (1/ms) cAMP hydrolysis rate by PDE
k_cAMP_PDEp     = 2*k_cAMP_PDE;   % (1/ms) cAMP hydrolysis rate by phosphorylated PDE
Km_PDE_cAMP     = 1.3;            % (uM) PDE affinity for cAMP

Kd_PDE_IBMX     = 30.0;           % (uM) Kd_R2cAMP_C for IBMX binding to PDE
k_PKA_PDE       = 7.5e-3;         % (1/ms) rate constant for PDE phosphorylation by type 1 PKA
k_PP_PDE        = 1.5e-3;         % (1/ms) rate constant for PDE dephosphorylation by phosphatases

% PKA module
PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA
PKItot          = 0.18;           % (uM) total PKI
kf_RC_cAMP      = 1;              % (1/[uM ms]) Kd for PKA RC binding to cAMP
kf_RCcAMP_cAMP  = 1;              % (1/[uM ms]) Kd for PKA RC:cAMP binding to cAMP
kf_RcAMPcAMP_C  = 4.375;          % (1/[uM ms]) Kd for PKA R:cAMPcAMP binding to C
kf_PKA_PKI      = 1;              % (1/[uM ms]) Ki for PKA inhibition by PKI
kr_RC_cAMP      = 1.64;           % (1/ms) Kd for PKA RC binding to cAMP
kr_RCcAMP_cAMP  = 9.14;           % (1/ms) Kd for PKA RC:cAMP binding to cAMP
kr_RcAMPcAMP_C  = 1;              % (1/ms) Kd for PKA R:cAMPcAMP binding to C
kr_PKA_PKI      = 2e-4;           % (1/ms) Ki for PKA inhibition by PKI
epsilon         = 10;             % (-) AKAP-mediated scaling factor

% PP1 module
PP1tot          = 0.89;           % (uM) total phosphatase 1
I1tot           = 0.3;            % (uM) total inhibitor 1
k_PKA_I1        = 60e-3;          % (1/ms) rate constant for I-1 phosphorylation by type 1 PKA
Km_PKA_I1       = 1.0;            % (uM) Km for I-1 phosphorylation by type 1 PKA
Vmax_PP2A_I1    = 14.0e-3;        % (uM/ms) Vmax for I-1 dephosphorylation by PP2A
Km_PP2A_I1      = 1.0;          	% (uM) Km for I-1 dephosphorylation by PP2A

Ki_PP1_I1       = 1.0e-3;         % (uM) Ki for PP1 inhibition by I-1
kf_PP1_I1       = 1;              % (uM) Ki for PP1 inhibition by I-1
kr_PP1_I1       = Ki_PP1_I1;      % (uM) Ki for PP1 inhibition by I-1

% LCC module
LCCtot          = 0.025;          % (uM) total ISO-type Ca channels
PKACII_LCCtot   = 0.025;          % (uM) type 2 PKA available to phosphorylate LCC
PP1_LCC         = 0.025;          % (uM) PP1 available to dephosphorylate LCC
PP2A_LCC        = 0.025;          % (uM) PP2A available to dephosphorylate LCC
k_PKA_LCC       = 54e-3;          % (1/ms) rate constant for LCC phosphorylation by type 2 PKA
Km_PKA_LCC      = 21;             % (uM) Km for LCC phosphorylation by type 2 PKA
k_PP1_LCC       = 8.52e-3;        % (1/ms) rate constant for LCC dephosphorylation by PP1
Km_PP1_LCC      = 3;            	% (uM) Km for LCC dephosphorylation by PP1
k_PP2A_LCC      = 10.1e-3;      	% (1/ms) rate constant for LCC dephosphorylation by PP2A
Km_PP2A_LCC     = 3;              % (uM) Km for LCC dephosphorylation by PP2A

% PLB module
PLBtot          = 106;            % (uM) total phospholamban
k_PKA_PLB       = 54e-3;          % (1/ms) rate constant for PLB phosphorylation by type 1 PKA
Km_PKA_PLB      = 21;             % (uM) Km for PLB phosphorylation by type 1 PKA
k_PP1_PLB       = 8.5e-3;         % (1/ms) rate constant for PLB dephosphorylation by PP1
Km_PP1_PLB      = 7.0;            % (uM) Km for PLB dephosphorylation by PP1

% PLM module
PLMtot          = 48;             % (uM) total phospholemman
k_PKA_PLM       = 54e-3;          % (1/ms) rate constant for PLM phosphorylation by type 1 PKA
Km_PKA_PLM      = 21;             % (uM) Km for PLM phosphorylation by type 1 PKA
k_PP1_PLM       = 8.5e-3;         % (1/ms) rate constant for PLM dephosphorylation by PP1
Km_PP1_PLM      = 7;              % (uM) Km for PLM dephosphorylation by PP1

% TnI module
TnItot          = 70;             % (uM) total troponin I
PP2A_TnI        = 0.67;           % (uM) PP2A available to dephosphorylate TnI
k_PKA_TnI       = 54e-3;          % (1/ms) rate constant for TnI phosphorylation by type 1 PKA
Km_PKA_TnI      = 21;             % (uM) Km for TnI phosphorylation by type 1 PKA
k_PP2A_TnI      = 10.1e-3;        % (1/ms) rate constant for TnI dephosphorylation by PP2A
Km_PP2A_TnI     = 4.1;            % (uM) Km for TnI dephosphorylation by PP2A

% ---- EC Coupling model parameters ------

% Universal parameters

R = 8314;   % R     [J/kmol*K]
Frdy = 96485;  % Frdy     [C/mol]

Vmyo = 20.8e-6;     % Vmyo  [uL]
Vnsr = 9.88e-7;     % Vnsr  [uL]
Vjsr = 9.3e-8;      % Vjsr  [uL]
ACap = 1.534e-4;    % ACap  [cm^2] with C = 1 uF/cm^2
Temp = 310;         % Temp  [K]
% extracellular concentrations     
Nao = 140;     % Extracellular Na  [mM]
Ko = 5.4;     % Extracellular K   [mM]
Cao = 1.8;     % Extracellular Ca  [mM]
% current conductances
g_Na = 8.0;    % G_Na      [mS/uF]
g_to = 0.35;   % G_to      [mS/uF] 
g_ss = 0.07;   % G_ss      [mS/uF]
g_kibar = 0.24;   % G_kibar   [mS/uF] 
g_kp = 0.008;  % G_kp      [mS/uF]
% I_Ca parameters
f = 300;    % f         [1/sec] 
g = 2e3;    % g         [1/sec]
gammao = 5187.5; % gammao    [1/sec/mM]
omega = 10;     % omega     [1/sec]
pCa = 5.823e-9*3.0;  % pCa       [cm/sec]
pK = 1.078e-11*3.0;  % pK        [cm/sec]
Nlcc = 3e5;    % Nlcc      [#/cell]
I_Ca05 = -0.458; % I_Ca05    [uA/uF]
% pumps and background currents
k_NaCa = 1483;   % k_NaCa    [uA/uF]
Km_Na = 87.5;   % Km_Na     [mM]
Km_Ca = 1.38;   % Km_Ca     [mM]
k_sat = 0.1;    % k_sat     [none]  
eta = 0.35;   % eta       [none]
ibarnak = 1.1;  % ibarnak   [uA/uF]
Km_Nai= 10;     % Km_Nai    [mM]
Km_Ko = 1.5;    % Km_Ko     [mM]
ibarpca = 1.15;   % ibarpca   [uA/uF]
Km_pca = 0.5e-3; % Km_pca    [mM]
g_Cab = 2.8e-3; % G_Cab     [uA/uF] 
g_Nab = 1.18e-3;   % G_Nab     [uA/uF] 
Pns = 0;      % Pns 
Km_ns = 1.2e-3; % Km_ns     [mM]
% Calcium handling parameters
I_upbar = 4.7;    % I_upbar   [mM/sec]
Km_up = 3e-4;   % Km_up     [mM]
nsrbar = 15;     % nsrbar    [mM]
tauon = 2e-3;   % tauon     [sec]
tauoff = 2e-3;   % tauoff    [sec]
gmaxrel = 60e3;   % gmaxrel   [mM/sec]
dcaith = 0.18e-3;% dcaith    [mM]
Km_rel = 0.8e-3; % Km_rel    [mM]
CSQNth = 8.75;   % CSQNth    [mM]
CSQNbar = 15;     % CSQNbar   [mM]
Km_csqn = 0.8;    % Km_csqn   [mM]
tau_tr = 5.7e-4; % tau_tr    [sec]
TRPNbar = 0.07;   % TRPNbar   [mM]
CMDNbar = 0.05;   % CMDNbar   [mM]
INDObar = 0.07;   % INDObar   [mM]
Km_trpn = 0.5128e-3;  % Km_trpn   [mM]
Km_cmdn = 2.38e-3;    % Km_cmdn   [mM]
Km_indo = 8.44e-4;    % Km_indo   [mM]
scaling_factor = factor;


%% Assemble Parameters Array
  params = [Ltot Atot FSK IBMX ...
   kf_Ra kr_Ra kf_LRa1 kr_LRa1 kf_LRi kr_LRi kf_LRa2 kr_LRa2...
  kf_RaG kr_RaG kf_LRaG1 kr_LRaG1 kf_LRaG2 kr_LRaG2 ...
  kf_ARa1 kr_ARa1 kf_ARi kr_ARi kf_ARa2 kr_ARa2 ... 
   kf_ARaG1 kr_ARaG1 kf_ARaG2 kr_ARaG2 ...
  b1ARtot Gstot kf_bARK kr_bARK k2_bARK Km_bARK kf_PKA kr_PKA k_G_act k_G_hyd...
  k_G_reassoc...
  ACtot ATP k_AC_basal Km_AC_basal ...
  kf_AC_Gsa kr_AC_Gsa k_AC_Gsa Km_AC_Gsa ...
  Kd_AC_FSK k_AC_FSK Km_AC_FSK ...
  PDEtot k_cAMP_PDE k_cAMP_PDEp Km_PDE_cAMP ...
  Kd_PDE_IBMX k_PKA_PDE k_PP_PDE ...
  PKAIItot PKItot ...
  kf_RC_cAMP kf_RCcAMP_cAMP kf_RcAMPcAMP_C kf_PKA_PKI ...
  kr_RC_cAMP kr_RCcAMP_cAMP kr_RcAMPcAMP_C kr_PKA_PKI ...
  epsilon ...
  PP1tot I1tot k_PKA_I1 Km_PKA_I1 Vmax_PP2A_I1 Km_PP2A_I1 ...
  kf_PP1_I1 kr_PP1_I1 ...
  LCCtot PKACII_LCCtot PP1_LCC PP2A_LCC ...
  k_PKA_LCC Km_PKA_LCC k_PP1_LCC Km_PP1_LCC k_PP2A_LCC Km_PP2A_LCC ...
  PLBtot k_PKA_PLB Km_PKA_PLB k_PP1_PLB Km_PP1_PLB ...
  PLMtot k_PKA_PLM Km_PKA_PLM k_PP1_PLM Km_PP1_PLM ...
  TnItot PP2A_TnI k_PKA_TnI Km_PKA_TnI k_PP2A_TnI Km_PP2A_TnI ...  
  Vmyo  Vnsr  Vjsr  ACap  Temp  Nao  Ko  Cao  g_Na  g_to  g_ss  ... 
  g_kibar  g_kp  f  g  gammao  omega  pCa  pK  Nlcc  I_Ca05  k_NaCa  Km_Na  Km_Ca  ...
  k_sat  eta  ibarnak  Km_Nai  Km_Ko  ibarpca  Km_pca  g_Cab  g_Nab  Pns  Km_ns  ...
  I_upbar  Km_up  nsrbar  tauon  tauoff  gmaxrel  dcaith  Km_rel  CSQNth  CSQNbar  ...
  Km_csqn  tau_tr  TRPNbar  CMDNbar  INDObar  Km_trpn  Km_cmdn  Km_indo...
  scaling_factor];

