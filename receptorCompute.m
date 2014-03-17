% Run file for receptor module of rat signaling model
%   RA: Robert Amanfu   <rka2p@virginia.edu>
%
% Robert Amanfu
% 11/08/11

%%
clear all; clc;
load KLcalc; load fitalpha;
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(18);
K_mod = KA;
KR = 10;
alpha_A = fitalpha(18);
gamma_A = 1;
KG = 2.4131;gamma_L =  0.3762;

PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA

RelTol = 1e-13;
MaxStep = 1e3;
options = odeset('MaxStep',MaxStep,'NonNegative',[1:2],'RelTol',RelTol);

params = receptorPARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,K_mod);
%% receptor species

load KLcalc; load fitalpha;
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(17);
alpha_A = fitalpha(17);
KG = 2.4131;gamma_L =  0.3762;
p = receptorPARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,K_mod);
CAR_dose = 0.3;
ISO_dose = 0;
y0 = zeros(2,1);[~,y] = ode15s(@receptorODE,[0; 20*60*1000],y0,options,p);
y0 = y(end,:);
tspan1 = [0;2*60*1000]; 
p(1) = 0;p(2) = 0;
[t1,y1] = ode15s(@receptorODE,tspan1,y0,options,p);
for tstep=1:length(t1),
    [~,algvars1(tstep,:)]=receptorODE(t1(tstep),y1(tstep,:),p);
end

tspan = [2*60*1000; 6*60*1000]; y0 = y1(end,:);
p(1) = 0.1;p(2) = CAR_dose;
[t2,y2] = ode15s(@receptorODE,tspan,y0,options,p);
for tstep=1:length(t2),
    [~,algvars2(tstep,:)]=receptorODE(t2(tstep),y2(tstep,:),p);
end

tspan = [6*60*1000; 10*60*1000]; y0 = y2(end,:);
p(1) = 10;p(2) = CAR_dose;
[t3,y3] = ode15s(@receptorODE,tspan,y0,options,p);
for tstep=1:length(t3),
    [~,algvars3(tstep,:)]=receptorODE(t3(tstep),y3(tstep,:),p);
end

tGly = [t1;t2;t3;]./60e3;
 yGly = [y1;y2;y3;];
algvarsGly = [algvars1;algvars2;algvars3;];
figure(1);

yCell=mat2cell(yGly,size(yGly,1),ones(size(yGly,2),1));
%%
[Ri,G] = yCell{:};
 Rtot = sum(algvarsGly,2) + Ri ;
% algvarsGly =  100*algvarsGly./repmat(Rtot,1,8);
 idx = Rtot< .0132+1e-4;
 algvarsCell=mat2cell(algvarsGly,size(algvarsGly,1),ones(size(algvarsGly,2),1));
[Ra, LRi ,LRa, RaG, LRaG ,ARi, ARa ,ARaG, ARii,ARai] =  algvarsCell{:};
%%
color2= [0.6 0 0];color1 = [0.5 0.5 0.5];
r = 0.0132;
 subplot(2,3,1);plot(tGly(idx),Ra(idx)./r,'LineWidth',2,'Color',color1);ylabel('Ra (% \beta1-AR)');hold all;xlabel('time (min)');
 subplot(2,3,2);plot(tGly(idx),ARii(idx)./r,'LineWidth',2,'Color',color1);ylabel('ARii (% \beta1-AR)');hold all;xlabel('time (min)');
subplot(2,3,3);plot(tGly(idx),RaG(idx)./r,'LineWidth',2,'Color',color1);ylabel('RaG (% \beta1-AR)');hold all;xlabel('time (min)');
subplot(2,3,4);plot(tGly(idx),LRaG(idx)./r,'LineWidth',2,'Color',color1);ylabel('LRaG (% \beta1-AR)');hold all;xlabel('time (min)');
subplot(2,3,5);plot(tGly(idx),ARi(idx)./r,'LineWidth',2,'Color',color1);ylabel('ARi (% \beta1-AR)');hold all;xlabel('time (min)');
subplot(2,3,6);plot(tGly(idx),Rtot(idx),'LineWidth',2,'Color',color1);xlabel('time (min)');ylabel('Rtot (\muM)');hold all;
%%
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperUnits', 'inches');
set(gcf, 'PaperPosition', [0.0 3.5 9.0  6.2]);
print -dpdf model_changes
