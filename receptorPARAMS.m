function params = receptorPARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,K_mod)

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


% ETC receptor module/Gs module
b1ARtot = 0.0132;
Gstot =  3.83;
kr_Ra  = KR;   
kr_LRa1  = (alpha_L*KL); 
kr_LRi = KL; 
kr_RaG = KG;
kr_LRaG2 = (gamma_L*KG);
kr_ARa1  = (alpha_A*KA);
kr_ARa1i = (alpha_A*K_mod);
kr_ARi = KA;
kr_ARii = K_mod;
kr_ARaG2 = (gamma_A*KG);





%% Assemble Parameters Array
  params =  [Ltot  Atot kr_Ra kr_LRa1 kr_LRi kr_RaG kr_LRaG2 kr_ARa1 ...
      kr_ARi  kr_ARaG2 b1ARtot Gstot kr_ARii kr_ARa1i];

