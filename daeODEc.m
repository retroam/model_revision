function [dydt,algvars] = daeODEc(t,y,params,flag)

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
  TnItot, PP2A_TnI, k_PKA_TnI, Km_PKA_TnI, k_PP2A_TnI, Km_PP2A_TnI, ...
 Vmyo, Vnsr, Vjsr, ACap, Temp, Nao, Ko ,Cao ,g_Na, g_to, g_ss, ... 
  g_kibar, g_kp, f, g, gammao, omega, pCa, pK, Nlcc, I_Ca05, k_NaCa, Km_Na, Km_Ca, ...
  k_sat, eta, ibarnak, Km_Nai, Km_Ko, ibarpca, Km_pca, g_Cab, g_Nab, Pns, Km_ns, ...
  I_upbar, Km_up, nsrbar, tauon, tauoff, gmaxrel, dcaith, Km_rel ,CSQNth, CSQNbar, ...
  Km_csqn, tau_tr, TRPNbar, CMDNbar, INDObar, Km_trpn, Km_cmdn, Km_indo,scaling_factor ] = pCell{:};

yCell = num2cell(y);
[Ri,G, b1AR_S464,b1AR_S301,...
    GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
  RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
  RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
  I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp, ...
    m, h, jo,v, w, x, yo, z,rto, sto ,ssto,rss ,...
    sss,Ca_nsr ,Ca_jsr,Nai ,Ki ,Cai,Vm,trelo]=yCell{:};



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

if flag == 1
    G_act = k_G_act*(RaG+LRaG + ARaG)*(1 - Atot/(scaling_factor + Atot));
else
    G_act = k_G_act*(RaG+LRaG + ARaG);
end
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
fracLCCap = LCCap/LCCtot;
fracLCCapo = 0.2041;


LCCb = LCCtot - LCCbp;
LCCb_phosph = epsilon*k_PKA_LCC*PKACII_LCC*LCCb/(Km_PKA_LCC+epsilon*LCCb);
LCCb_dephosph = epsilon*k_PP1_LCC*PP1_LCC*LCCbp/(Km_PP1_LCC+epsilon*LCCbp);
dLCCbp = LCCb_phosph - LCCb_dephosph;
fracLCCbp = LCCbp/LCCtot;
fracLCCbpo = 0.2336;

% PLB module
PLB = PLBtot - PLBp;
PLB_phosph = k_PKA_PLB*PKACI*PLB/(Km_PKA_PLB+PLB);
PLB_dephosph = k_PP1_PLB*PP1*PLBp/(Km_PP1_PLB+PLBp);
dPLBp = PLB_phosph - PLB_dephosph;
fracPLBp = PLBp/PLBtot;
fracPLB = PLB/PLBtot;
fracPLBo = 0.9613; % adjust when changes are made to signaling model!

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

% -------- EC COUPLING MODEL -----------

% Constants
R = 8314;   % R     [J/kmol*K]
Frdy = 96485;  % Frdy     [C/mol]
FoRT = Frdy/R/Temp;
zna = 1;    % Na valence
zk = 1;     % K valence
zca = 2;    % Ca valence
% Nernst Potentials
ena = (1/FoRT/zna)*log(Nao/Nai);       % should be 70.54 mV
ek = (1/FoRT/zk)*log(Ko/Ki);		 % should be -87.94 mV
eca = (1/FoRT/zca)*log(Cao/Cai);  % should be 120 mV
%eks = (1/FoRT)*logn((ko+prnak*nao)/(ki+prnak*nai))	% should be -77.54
ecl = -40.0;

% I_Na: Fast Na Current
am = 0.32*(Vm+47.13)/(1-exp(-0.1*(Vm+47.13)));
bm = 0.08*exp(-Vm/11);
if Vm >= -40
    ah = 0; aj = 0;
    bh = 1/(0.13*(1+exp(-(Vm+10.66)/11.1)));
    bj = 0.3*exp(-2.535e-7*Vm)/(1+exp(-0.1*(Vm+32)));
else
    ah = 0.135*exp((80+Vm)/-6.8);
    bh = 3.56*exp(0.079*Vm)+3.1e5*exp(0.35*Vm);
    aj = (-1.2714e5*exp(0.2444*Vm)-3.474e-5*exp(-0.04391*Vm))*(Vm+37.78)/(1+exp(0.311*(Vm+79.23)));
    bj = 0.1212*exp(-0.01052*Vm)/(1+exp(-0.1378*(Vm+40.14)));
end
dm = 1e3*(am*(1-m)-bm*m);
dh = 1e3*(ah*(1-h)-bh*h);
djo = 1e3*(aj*(1-jo)-bj*jo);
I_Na = g_Na*m^3*h*jo*(Vm-ena);

% I_Ca: L-type Calcium Current
alcc = 400*exp( (Vm+2)/10 );    % [1/sec]
blcc = 50*exp( -1*(Vm+2)/13 );  % [1/sec]
flcc = f*(0.375*fracLCCap/fracLCCapo+0.625); % PHOSPHOREGULATION
ylccinf = 1.0/(1+exp((Vm+55.0)/7.5 )) +0.1/(1+exp((-Vm+21.0)/6.0 ));
tauylcc = 0.02 + 0.3/( 1 +exp( (Vm+30)/9.5 ) );    % [sec]			
gamma = gammao*Cai; % [1/sec/mM]		
vgamma = gamma*((1-v)^4+2*v*(1-v)^3+4*v^2*(1-v)^2+8*v^3*(1-v)+16*v^4*(1-flcc/g ));
vomega = omega*((1-w)^4+1/2*w*(1-w)^3+1/4*w^2*(1-w)^2+1/8*w^3*(1-w)+1/16*w^4);

dv = alcc*(1-v)-blcc*v;
dw = 2*alcc*(1-w)-blcc/2*w;
dx = flcc*(1-x)-g*x;
dyo = (ylccinf - yo)/tauylcc;
dz = vomega*(1-z)-vgamma*z;

if abs(Vm) < 1e-6
    disp('Warning! Voltage near zero, could influence ibarca and ibark!');
end
ibarca = pCa*4*(Vm*Frdy*FoRT) * (1e-3*exp(2*Vm*FoRT)-0.341*Cao) /(exp(2*Vm*FoRT)-1);
ibark = pK*(Vm*Frdy*FoRT)*(Ki*exp(Vm*FoRT)-Ko) /(exp(Vm*FoRT)-1);
favail = 0.5*(0.4*fracLCCbp/fracLCCbpo+0.60);   % PHOSPHOREGULATION
I_Ca = ibarca*Nlcc *favail*v^4*x*yo*z;
I_CaK = ibark/(1+I_Ca/I_Ca05 )*Nlcc *favail*v^4*x*yo*z;
%I_Ca = 0; I_CaK = 0;
I_Catot = I_Ca+I_CaK;


% I_to: Transient Outward K Current
rtoss = 1/(1+exp((Vm+10.6)/-11.42));
stoss = 1/(1+exp((Vm+45.3)/6.8841));
taurto = 1.0/(45.16*exp(0.03577*(Vm+50))+98.9*exp(-0.1*(Vm+38))); % [sec]
tausto = 0.35*exp(-((Vm+70)/15)^2)+0.035;	    % [sec]		 
taussto = 3.7*exp(-((Vm+70)/30)^2)+0.035;        % [sec]			 
drto = (rtoss-rto)/taurto;
dsto = (stoss-sto)/tausto;
dssto = (stoss-ssto)/taussto;
I_to = g_to*rto*(0.886*sto+0.114*ssto)*(Vm-ek);   % [uA/uF]

% I_ss: Steady-state K Current
rssinf = 1/(1+exp(-(Vm+11.5)/11.82));
taurss = 10/(45.16*exp(0.03577*(Vm+50))+98.9*exp(-0.1*(Vm+38)));
sssinf = 1/(1+exp((Vm+87.5)/10.3));
tausss = 2.1;
drss = (rssinf-rss)/taurss;
dsss = (sssinf-sss)/tausss;
I_ss = g_ss*rss*sss*(Vm-ek);

% I_ki: Time-Independent K Current
aki = 1.02/(1+exp(0.2385*(Vm-ek-59.215)));
bki =(0.49124*exp(0.08032*(Vm+5.476-ek)) + exp(0.06175*(Vm-ek-594.31))) /(1 + exp(-0.5143*(Vm-ek+4.753)));
kiss = aki/(aki+bki);
I_ki = g_kibar*sqrt(Ko/5.4)*kiss*(Vm-ek) ;

% I_kp: Plateau K Current
kp = 1/(1+exp((7.488-Vm)/5.98));
I_kp = g_kp*kp*(Vm-ek);

% I_ncx: Na/Ca Exchanger Current
s4 = exp(eta*Vm*FoRT)*Nai^3*Cao;
s5 = exp((eta-1)*Vm*FoRT)*Nao^3*Cai;
I_ncx = k_NaCa/(Km_Na^3+Nao^3) /(Km_Ca+Cao) /(1+k_sat*exp((eta-1)*Vm*FoRT)) *(s4-s5);

% I_nak: Na/K Pump Current
sigma = (exp(Nao/67.3)-1)/7;
fnak = 1/(1+0.1245*exp(-0.1*Vm*FoRT)+0.0365*sigma*exp(-Vm*FoRT));
I_nak = ibarnak *fnak /(1+(Km_Nai/Nai)^1.5) *(Ko/(Ko+Km_Ko));

% I_pca: Sarcolemmal Ca Pump Current
I_pca = ibarpca*Cai/(Km_pca+Cai);
 
% I_cab: Ca Background Current
I_cab = g_Cab*(Vm-eca);
 
% I_nab: Na Background Current
I_nab = g_Nab*(Vm-ena);

% I_nsca: Nonspecific Ca-Activated Current: not used

% Total Membrane Currents
I_Na_tot = I_Na+I_nab+3*I_ncx+3*I_nak;          % [uA/uF]
I_K_tot = I_to+I_ss+I_ki+I_kp-2*I_nak+I_CaK;    % [uA/uF]
I_Ca_tot = I_Ca+I_cab+I_pca-2*I_ncx;            % [uA/uF]

% Calcium Induced Calcium Release (CICR)
trel = trelo+2e-3;
ryron = 1-exp(-trel/tauon);  
ryroff = exp(-trel/tauoff);
grel = gmaxrel/(1+exp((I_Ca_tot+5)/0.9));         % adjusted for rat  [1/sec] 
I_rel = grel*ryron*ryroff*(Ca_jsr -Cai);        %   [mM/sec]

% Other SR fluxes and concentrations
Km_up = Km_up*(1+2*fracPLB)/(1+2*fracPLBo);     % PHOSPHOREGULATION
I_up = I_upbar*Cai^2/(Km_up^2+Cai^2);         %   [mM/sec]
% Original: I_up = I_upbar*Cai^2/(Km_up^2+Cai^2);         %   [mM/sec]
I_leak = I_upbar*Ca_nsr /nsrbar ;                     %   [mM/sec]
I_tr = (Ca_nsr -Ca_jsr )/tau_tr ;                    %   [mM/sec]
Bjsr = 1/( 1+CSQNbar  *Km_csqn /(Km_csqn +Ca_jsr )^2 );
dCa_nsr  = I_up-I_leak-I_tr*Vjsr/Vnsr;        %   [mM/sec]
dCa_jsr  = Bjsr*(I_tr-I_rel);                   %   [mM/sec]
SRcontent = 1e3*((Ca_jsr +Ca_jsr /Bjsr)*Vjsr/Vmyo+Ca_nsr *Vnsr/Vmyo);    % [umol/L cytosol]

% Cytoplasmic Calcium Buffering
btrpn = TRPNbar*Km_trpn/(Km_trpn+Cai)^2;
bcmdn = CMDNbar*Km_cmdn /(Km_cmdn +Cai)^2;
bindo = INDObar*Km_indo /(Km_indo +Cai)^2;
Bmyo = 1/( 1+ bcmdn + btrpn + btrpn + bindo);

% Ion Concentrations and Membrane Potential
dNai = -1e3*I_Na_tot*ACap/(Vmyo*zna*Frdy);          % [mM/sec] 
dKi = -1e3*I_K_tot*ACap/(Vmyo*zk*Frdy);            % [mM/sec]
dCai = -Bmyo*(1e3*I_Ca_tot*ACap/(Vmyo*zca*Frdy) ...
    +((I_up-I_leak)*Vnsr/Vmyo)-(I_rel*Vjsr/Vmyo));    % [mM/sec]

% Simulation type
protocol = 'pace';

switch lower(protocol)
    case {'none',''},
        I_app = 0;
    case 'pace',        % pace w/ current injection at rate 'rate'
		rate = 1;
		if mod(t+0.9,1/rate) <= 5e-3
            I_app = 10.0;
		else
            I_app = 0.0;
		end
    case 'vclamp',      
		V_hold = -40;
        V_test = -10;
		if (t > 59.1 & t < 59.5)
		    V_clamp = V_test;
		else
		    V_clamp = V_hold;
		end
		R_clamp = .02; 
		I_app = (V_clamp-Vm)/R_clamp;
end  

dVm = -1e3*(I_Ca_tot+I_K_tot+I_Na_tot-I_app);

% CICR timing: trelo tracks the last time when Vdot > 30 mV/msec
if (dVm > 30e3) 
    dtrelo = 1-1e4*trelo; 
else
    dtrelo = 1;
end

% ----- END EC COUPLING MODEL ---------------





%% Reassemble dydt


dydt1 = [dRi dG db1AR_S464...
   db1AR_S301 dGsaGTPtot dGsaGDP dGsby ...
  dAC_GsaGTP dcAMPtot dPDEp ...
  dRC_I dRCcAMP_I dRCcAMPcAMP_I dRcAMPcAMP_I dPKACI dPKACI_PKI ...
  dRC_II dRCcAMP_II dRCcAMPcAMP_II dRcAMPcAMP_II dPKACII dPKACII_PKI ...
  dI1p_PP1 dI1ptot dLCCap dLCCbp dPLBp dPLMp dTnIp]';

 dydt2 = [dm dh djo dv dw dx dyo dz drto dsto dssto drss ...
     dsss dCa_nsr dCa_jsr dNai dKi dCai dVm dtrelo ]';

dydt = [dydt1*1e3 ;dydt2];
algvars = [Ra LRi LRa RaG LRaG ARi ARa ARaG];
