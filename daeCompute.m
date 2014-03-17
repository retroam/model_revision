% Run file for rat signaling model
% 
%
% Copyright 2011, Cardiac Systems Biology Lab, University of Virginia
%   JS: Jeff Saucerman  <jsaucerman@virginia.edu>
%   JY: Jason Yang      <jhyang@virginia.edu>
%   RA: Robert Amanfu   <rka2p@virginia.edu>
%
% Robert Amanfu
% 11/08/11
%% Receptor parameters to optimize
 clc; clear all;%close all;
           
 
PKAItot = 0.59;           % (uM) total type 1 PKA
PKAIItot = 0.025;          % (uM) total type 2 PKA
KR = 10;
gamma_A = 1;


%% Checking receptor species
% 
scaling_factor = 0.1;
load KLcalc; load fitalpha;
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(17);
alpha_A = fitalpha(17);
KG = 2.4131;gamma_L =  0.3762;
RelTol = 1e-13;
MaxStep = 1e3;
options = odeset('MaxStep',MaxStep,'NonNegative',[1:29],'RelTol',RelTol);
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,scaling_factor);

flag = 0;
y0 = zeros(29,1);[~,y] = ode15s(@daeODE,[0; 20*60*1000],y0,options,p,flag);
y0 = y(end,:);

tspan1 = [0;2*60*1000]; 
p(1) = 0;p(2) = 0;flag = 0;
[t1,y1] = ode15s(@daeODE,tspan1,y0,options,p,flag);
for tstep=1:length(t1),
    [~,algvars1(tstep,:)]=daeODE(t1(tstep),y1(tstep,:),p,flag);
end

tspan = [2*60*1000; 2.1*60*1000]; y0 = y1(end,:);
p(1) = 0.1;p(2) = 0.1;flag = 1;
[t2,y2] = ode15s(@daeODE,tspan,y0,options,p,flag);
for tstep=1:length(t2),
    [~,algvars2(tstep,:)]=daeODE(t2(tstep),y2(tstep,:),p,flag);
end


tspan = [2.1*60*1000; 6*60*1000]; y0 = y2(end,:);
p(1) = 0.1;p(2) = 0.1;flag = 1;
[t3,y3] = ode15s(@daeODE,tspan,y0,options,p, flag);

for tstep=1:length(t3),
    [~,algvars3(tstep,:)]=daeODE(t3(tstep),y3(tstep,:),p,flag);
end
tspan = [6*60*1000; 10*60*1000]; y0 = y3(end,:);
p(1) = 10;p(2) = 0.1;flag = 1;
[t4,y4] = ode15s(@daeODE,tspan,y0,options,p,flag);
for tstep=1:length(t4),
    [~,algvars4(tstep,:)]=daeODE(t4(tstep),y4(tstep,:),p,flag);
end

tGly = [t1;t2;t3;t4]./60e3;
 yGly = [y1;y2;y3;y4];
algvarsGly = [algvars1;algvars2;algvars3;algvars4];
yCell=mat2cell(yGly,size(yGly,1),ones(size(yGly,2),1));

[Ri,G, b1AR_S464,b1AR_S301,...
    GsaGTPtot, GsaGDP, Gsby, AC_GsaGTP ,cAMPtot, PDEp, ...
  RC_I, RCcAMP_I, RCcAMPcAMP_I, RcAMPcAMP_I, PKACI ,PKACI_PKI, ...
  RC_II, RCcAMP_II, RCcAMPcAMP_II, RcAMPcAMP_II, PKACII, PKACII_PKI, ...
  I1p_PP1, I1ptot, LCCap ,LCCbp, PLBp, PLMp, TnIp] = yCell{:};
 Rtot = sum(algvarsGly,2) + b1AR_S464 + b1AR_S301 ;
 index = Rtot> .0132+1e-4;algvarsGly(index,:) = NaN;talgvars = tGly; talgvars(index) = NaN;
 algvarsCell=mat2cell(algvarsGly,size(algvarsGly,1),ones(size(algvarsGly,2),1));
[Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG] =  algvarsCell{:};

subplot(3,2,1);plot(talgvars(~isnan(talgvars)),Ra(~isnan(Ra)));title('Ra');axis tight;hold all;
subplot(3,2,2);plot(talgvars(~isnan(talgvars)),LRa(~isnan(LRa)));title('LRa');axis tight;hold all;
subplot(3,2,3);plot(talgvars(~isnan(talgvars)),RaG(~isnan(RaG)));title('RaG');axis tight;hold all;
subplot(3,2,4);plot(talgvars(~isnan(talgvars)),LRaG(~isnan(LRaG)));title('LRaG');axis tight;hold all;
subplot(3,2,5);plot(talgvars(~isnan(talgvars)),ARi(~isnan(LRaG)));title('ARi');axis tight;hold all;
subplot(3,2,6);plot(tGly,cAMPtot);title('cAMPtot');axis tight;hold all;
xlabel('time (mins)');ylabel('TnIp(\muM)');hold all;