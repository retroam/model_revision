function params = receptorPARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A)

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
  ];

