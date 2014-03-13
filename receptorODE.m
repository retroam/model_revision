function [dydt,algvars] = receptorODE(t,y,params)

% ODE for rat b-adrenergic signaling model(based on Jason Yang's mouse
% model)
% dydt = ratODE(t,y,params)
%
% Copyright 2011, Cardiac Systems Biology Lab, University of Virginia
%   JS: Jeff Saucerman  <jsaucerman@virginia.edu>
%   JY: Jason Yang      <jhyang@virginia.edu>
%   RA: Robert Amanfu   <rka2p@virginia.edu>
%
% Robert Amanfu
% 11/08/11

%% Assign Parameters and State Variables
pCell = num2cell(params);
  [Ltot ,Atot, FSK, IBMX, ...
   kf_Ra, kr_Ra, kf_LRa1, kr_LRa1, kf_LRi,kr_LRi ,kf_LRa2, kr_LRa2,...
  kf_RaG ,kr_RaG, kf_LRaG1, kr_LRaG1, kf_LRaG2, kr_LRaG2, ...
  kf_ARa1, kr_ARa1, kf_ARi, kr_ARi, kf_ARa2, kr_ARa2, ... 
   kf_ARaG1 ,kr_ARaG1, kf_ARaG2, kr_ARaG2,  ...
  b1ARtot, Gstot ,kf_bARK, kr_bARK,k2_bARK ,Km_bARK, kf_PKA ,kr_PKA, k_G_act, k_G_hyd,...
  k_G_reassoc,...
  ACtot ,ATP, k_AC_basal ,Km_AC_basal, ...
  kf_AC_Gsa, kr_AC_Gsa, k_AC_Gsa, Km_AC_Gsa, ...
  Kd_AC_FSK, k_AC_FSK, Km_AC_FSK ,...
  PDEtot, k_cAMP_PDE ,k_cAMP_PDEp, Km_PDE_cAMP, ...
  Kd_PDE_IBMX, k_PKA_PDE ,k_PP_PDE, ...
  PKAIItot, PKItot, ...
  kf_RC_cAMP, kf_RCcAMP_cAMP, kf_RcAMPcAMP_C, kf_PKA_PKI, ...
  kr_RC_cAMP, kr_RCcAMP_cAMP, kr_RcAMPcAMP_C, kr_PKA_PKI, ...
  epsilon, ...
  PP1tot, I1tot, k_PKA_I1, Km_PKA_I1, Vmax_PP2A_I1, Km_PP2A_I1, ...
  kf_PP1_I1, kr_PP1_I1, ...
  LCCtot, PKACII_LCCtot ,PP1_LCC, PP2A_LCC, ...
  k_PKA_LCC, Km_PKA_LCC, k_PP1_LCC ,Km_PP1_LCC, k_PP2A_LCC, Km_PP2A_LCC,...
  PLBtot, k_PKA_PLB, Km_PKA_PLB, k_PP1_PLB, Km_PP1_PLB, ...
  PLMtot, k_PKA_PLM ,Km_PKA_PLM, k_PP1_PLM, Km_PP1_PLM, ...
  TnItot, PP2A_TnI, k_PKA_TnI, Km_PKA_TnI, k_PP2A_TnI, Km_PP2A_TnI ...
  ] = pCell{:};

yCell = num2cell(y);
[Ri,G, b1AR_S464,b1AR_S301,...
    GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
  RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
  RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
  I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp]=yCell{:};



%% Signaling Module

% b-AR/Gs module

% ETC module
 
Ra = Ri/kr_Ra; % extended ternary complex model
LRi = Ltot*Ri/kr_LRi;
LRa = Ltot*Ra/kr_LRa1;
RaG = Ra*G/kr_RaG;
LRaG = LRa*G/kr_LRaG2;
ARi = Atot*Ri/kr_ARi;
ARa = Atot*Ra/kr_ARa1;
ARaG = ARa*G/kr_ARaG2;
b1ARact = b1ARtot - b1AR_S464 - b1AR_S301;
dRi = b1ARact -  Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG - Ri;
dG = Gstot - LRaG - RaG - ARaG  - G;

bARK_desens =  kf_bARK*(LRa+LRaG+RaG+ ARa + ARaG);
bARK_resens = kr_bARK*b1AR_S464;
PKA_desens =  kf_PKA*PKACI*b1ARact;
PKA_resens =  kr_PKA*b1AR_S301;
db1AR_S464 = bARK_desens - bARK_resens;
db1AR_S301 = PKA_desens - PKA_resens;

G_act = k_G_act*(RaG+LRaG + ARaG);
G_hyd = k_G_hyd*GsaGTPtot;
G_reassoc = k_G_reassoc*GsaGDP*Gsby;
dGsaGTPtot = G_act - G_hyd;
dGsaGDP = G_hyd - G_reassoc;
dGsby = G_act - G_reassoc;

% cAMP module
cAMP = cAMPtot - (RCcAMP_I+2*RCcAMPcAMP_I+2*RcAMPcAMP_I) - (RCcAMP_II+2*RCcAMPcAMP_II+2*RcAMPcAMP_II);
AC = ACtot-AC_GsaGTP;
GsaGTP = GsaGTPtot - AC_GsaGTP;
dAC_GsaGTP = kf_AC_Gsa*GsaGTP*AC - kr_AC_Gsa*AC_GsaGTP;

AC_FSK = FSK*AC/Kd_AC_FSK;
AC_ACT_BASAL = 0;
AC_ACT_GSA = k_AC_Gsa*AC_GsaGTP*ATP/(Km_AC_Gsa+ATP);
AC_ACT_FSK = k_AC_FSK*AC_FSK*ATP/(Km_AC_FSK+ATP);

PDE_IBMX = PDEtot*IBMX/Kd_PDE_IBMX;
PDE = PDEtot - PDE_IBMX - PDEp;
dPDEp = k_PKA_PDE*PKACII*PDE - k_PP_PDE*PDEp;
PDE_ACT = k_cAMP_PDE*PDE*cAMP/(Km_PDE_cAMP+cAMP) + k_cAMP_PDEp*PDEp*cAMP/(Km_PDE_cAMP+cAMP);

dcAMPtot = AC_ACT_BASAL + AC_ACT_GSA + AC_ACT_FSK - PDE_ACT;


% PKA module
PKI = PKItot - PKACI_PKI - PKACII_PKI;

dRC_I = - kf_RC_cAMP*RC_I*cAMP + kr_RC_cAMP*RCcAMP_I;
dRCcAMP_I = - kr_RC_cAMP*RCcAMP_I + kf_RC_cAMP*RC_I*cAMP - kf_RCcAMP_cAMP*RCcAMP_I*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_I;
dRCcAMPcAMP_I = - kr_RCcAMP_cAMP*RCcAMPcAMP_I + kf_RCcAMP_cAMP*RCcAMP_I*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_I + kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI;
dRcAMPcAMP_I = - kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I;
dPKACI = -kr_RcAMPcAMP_C*RcAMPcAMP_I*PKACI + kf_RcAMPcAMP_C*RCcAMPcAMP_I - kf_PKA_PKI*PKACI*PKI + kr_PKA_PKI*PKACI_PKI;
dPKACI_PKI = - kr_PKA_PKI*PKACI_PKI + kf_PKA_PKI*PKACI*PKI;

dRC_II = - kf_RC_cAMP*RC_II*cAMP + kr_RC_cAMP*RCcAMP_II;
dRCcAMP_II = - kr_RC_cAMP*RCcAMP_II + kf_RC_cAMP*RC_II*cAMP - kf_RCcAMP_cAMP*RCcAMP_II*cAMP + kr_RCcAMP_cAMP*RCcAMPcAMP_II;
dRCcAMPcAMP_II = - kr_RCcAMP_cAMP*RCcAMPcAMP_II + kf_RCcAMP_cAMP*RCcAMP_II*cAMP - kf_RcAMPcAMP_C*RCcAMPcAMP_II + kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII;
dRcAMPcAMP_II = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II;
dPKACII = - kr_RcAMPcAMP_C*RcAMPcAMP_II*PKACII + kf_RcAMPcAMP_C*RCcAMPcAMP_II - kf_PKA_PKI*PKACII*PKI + kr_PKA_PKI*PKACII_PKI;
dPKACII_PKI = - kr_PKA_PKI*PKACII_PKI + kf_PKA_PKI*PKACII*PKI;

% I-1/PP1 module
I1 = I1tot - I1ptot;
PP1 = PP1tot - I1p_PP1;
I1p = I1ptot - I1p_PP1;
I1_phosph = k_PKA_I1*PKACI*I1/(Km_PKA_I1+I1);
I1_dephosph = Vmax_PP2A_I1*I1ptot/(Km_PP2A_I1+I1ptot);

dI1p_PP1 = kf_PP1_I1*PP1*I1p - kr_PP1_I1*I1p_PP1;
dI1ptot = I1_phosph - I1_dephosph;

% LCC module
PKACII_LCC = (PKACII_LCCtot/PKAIItot)*PKACII;
LCCa = LCCtot - LCCap;
LCCa_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCa/(Km_PKA_LCC+epsilon*LCCa);
LCCa_dephosph = epsilon*k_PP2A_LCC*PP2A_LCC*LCCap/(Km_PP2A_LCC+epsilon*LCCap);
dLCCap = LCCa_phosph - LCCa_dephosph;

LCCb = LCCtot - LCCbp;
LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb);
LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp);
dLCCbp = LCCb_phosph - LCCb_dephosph;

% PLB module
PLB = PLBtot - PLBp;
PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB);
PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp);
dPLBp = PLB_phosph - PLB_dephosph;

% PLM module
PLM = PLMtot - PLMp;
PLM_phosph = k_PKA_PLM*PKACI*PLM/(Km_PKA_PLM+PLM);
PLM_dephosph = k_PP1_PLM*PP1*PLMp/(Km_PP1_PLM+PLMp);
dPLMp = PLM_phosph - PLM_dephosph;

% TnI module
TnI = TnItot - TnIp;
TnI_phosph = k_PKA_TnI*PKACI*TnI/(Km_PKA_TnI+TnI);
TnI_dephosph = k_PP2A_TnI*PP2A_TnI*TnIp/(Km_PP2A_TnI+TnIp);
dTnIp = TnI_phosph - TnI_dephosph;



%% Reassemble dydt


dydt = [dRi dG db1AR_S464...
   db1AR_S301 dGsaGTPtot dGsaGDP dGsby ...
  dAC_GsaGTP dcAMPtot dPDEp ...
  dRC_I dRCcAMP_I dRCcAMPcAMP_I dRcAMPcAMP_I dPKACI dPKACI_PKI ...
  dRC_II dRCcAMP_II dRCcAMPcAMP_II dRcAMPcAMP_II dPKACII dPKACII_PKI ...
  dI1p_PP1 dI1ptot ...
  dLCCap dLCCbp dPLBp dPLMp dTnIp ...
 ]';

algvars = [Ra LRi LRa RaG LRaG ARi ARa ARaG];
