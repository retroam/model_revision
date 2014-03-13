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
  [Ltot ,Atot,kr_Ra, kr_LRa1,kr_LRi,kr_RaG, kr_LRaG2, kr_ARa1,...
      kr_ARi, kr_ARaG2, b1ARtot, Gstot] = pCell{:};

yCell = num2cell(y);
[Ri,G]=yCell{:};



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
b1ARact = b1ARtot;
dRi = b1ARact -  Ra - LRi - LRa - RaG - LRaG - ARi - ARa - ARaG - Ri;
dG = Gstot - LRaG - RaG - ARaG  - G;





%% Reassemble dydt


dydt = [dRi dG]';

algvars = [Ra LRi LRa RaG LRaG ARi ARa ARaG];
