 
close all; clear all;clc;
KL = 0.2;
KR =  10;%10 [uM] 
alpha_L = 1/32;% 1/32
KA = 500e-6;%500e-6 2e-3;
alpha_A = 1; %1 1.5574
gamma_A =  1; %1 
PKAItot         = 0.59;           % (uM) total type 1 PKA
PKAIItot        = 0.025;          % (uM) total type 2 PKA


load fitalpha.mat;
 drugList = {'Isoproterenol'; 'Basal' ;'Epinephrine'; 'Norepinephrine'; 'Fenoterol' ;'Formoterol';
'Salbutamol' ;'Terbutaline' ;'Broxaterol' ;'Salmeterol'; 'BRL-37344'; 'CGP-12177';
'Alprenolol' ;'Pindolol'; 'SR 59230A'; 'Atenolol' ;'Carvedilol' ;'Metoprolol'; 'Bisoprolol';
'Propranolol' ;'CGP-20712' ;'ICI-118551'};

 Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
      5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
      1.80e-03,4.50e-03,4.95e-02];
  KLcalc =Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));
   KLcalc = (0.2/KLcalc(1))*KLcalc; 
KG = 2.4131;gamma_L =  0.3762;
flag = 0;
scaling_factor = .326;
%%
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,scaling_factor);

RelTol = 1e-13;MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p,flag);y0 = y(end,:);
tspan = [0;2*60*1000]; 
dose = 10.^[-5:.1:3];

 beta_blocker = [20 18 17];
 flag = [0 0 1];
% beta_blocker = [1:22];
for i=1:length(beta_blocker)
    for j = 1:length(dose)
        disp(['Beta Blocker:' drugList{beta_blocker(i)} ' dose: ' num2str(dose(j))]);
        KL =  KLcalc(1);
        alpha_L = fitalpha(1);
        KA = KLcalc(beta_blocker(i)); alpha_A = fitalpha(beta_blocker(i));
        p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,scaling_factor);
        tspan = [0;2*60*1000];
        p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
        [t1,y] = ode15s(@daeODE,tspan,y0,options,p,flag(1));
        bg(j,1) = max(y(end,9));
        yfinal = y(end,:); 
        tspan = [2*60*1000;6*60*1000];
        p(1) = 0.1; p(2) = dose(j);  % agonist conc and beta blocker conc
        [t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p,flag(i)); 
        bg(j,2) = max(y2(end,9));
        yfinal = y2(end,:); 
        tspan2 = [6*60*1000 10*60*1000]; %changed the time
        p(1) = 10; p(2) = dose(j);     % agonist conc and beta blocker conc
        [t3,y3]=ode15s(@daeODE,tspan2,yfinal,options,p,flag(i));
        bg(j,3) = max(y3(end,9));
    end 
         sensitivity(:,i) = (bg(:,3)- bg(:,2));
           response(:,i) = (bg(:,2));

end
%%
close all
figure;semilogx(dose*1e-6,sensitivity,'Linewidth',2);legend('propranolol','metoprolol','carvedilol')
xlabel('dose (M)');ylabel('cAMP sensitivity (\muM)');box off;xlim([1e-12 1e-2])
line('XData',[1e-7 1e-7],'YData',[-1 3],'LineStyle','-.','LineWidth',2);
line('XData',[1e-6 1e-6],'YData',[-1 3],'LineStyle','-.','LineWidth',2);
line('XData',[0.3e-6 0.3e-6],'YData',[-1 3],'LineStyle','-.','LineWidth',2);
set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.9 4.4 2.8 2.3]);
cd('A:\Robert\Manuscripts\Modeling Beta1-adrenergic receptor blockers and polymorphisms in cardiac myocytes\Figures\Reviewers');
print -dpdf BBdoseAwflag
%%
close all
figure;subplot(1,3,1);semilogx(dose*1e-6,[sensitivity(:,1) response(:,1)],'Linewidth',2);hold all;
xlabel('dose (M)');ylabel('cAMP(\muM)');box off;xlim([1e-12 1e-2]);title('PRO');
set(gca,'Xtick',[1e-010 ,1e-06,0.01])
line('XData',[1e-7 1e-7],'YData',[-1 4],'LineStyle','-.','LineWidth',2);
subplot(1,3,2);semilogx(dose*1e-6,[sensitivity(:,2) response(:,2)],'Linewidth',2);
xlabel('dose (M)');ylabel('cAMP(\muM)');box off;xlim([1e-12 1e-2]);title('MET');
set(gca,'Xtick',[1e-010 ,1e-06,0.01])
line('XData',[1e-6 1e-6],'YData',[-1 4],'LineStyle','-.','LineWidth',2);
subplot(1,3,3);semilogx(dose*1e-6,[sensitivity(:,3) response(:,3)],'Linewidth',2);
xlabel('dose (M)');ylabel('cAMP(\muM)');box off;xlim([1e-12 1e-2]);title('CAR');
line('XData',[0.3e-6 0.3e-6],'YData',[-1 4],'LineStyle','-.','LineWidth',2);
line('XData',[1e-6 1e-6],'YData',[-1 4],'LineStyle','-.','LineWidth',2);
 set(gca,'Xtick',[1e-010 ,1e-06,0.01])
 set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [0.5 4.4 7.4 2.3]);
cd('A:\Robert\Manuscripts\Modeling Beta1-adrenergic receptor blockers and polymorphisms in cardiac myocytes\Figures\Reviewers');
print -dpdf BBdoseBwflag

%% carvedilol with change in alpha

for i=1:length(beta_blocker)
    for j = 1:length(dose)
        disp(['Beta Blocker:' drugList{beta_blocker(i)} ' dose: ' num2str(dose(j))]);
        KL =  KLcalc(1);
        alpha_L = fitalpha(1);
        KA = KLcalc(beta_blocker(i)); alpha_A = fitalpha(beta_blocker(i));
        p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
        tspan = [0;2*60*1000];
        p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
        [t1,y] = ode15s(@daeODE,tspan,y0,options,p);
        bg(j,1) = max(y(end,9));
        yfinal = y(end,:); 
        tspan = [2*60*1000;6*60*1000];
        p(1) = 0.1; p(2) = dose(j);  % agonist conc and beta blocker conc
        [t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
        bg(j,2) = max(y2(end,9));
        yfinal = y2(end,:); 
        tspan2 = [6*60*1000 10*60*1000]; %changed the time
        p(1) = 10; p(2) = dose(j);     % agonist conc and beta blocker conc
        [t3,y3]=ode15s(@daeODE,tspan2,yfinal,options,p);
        bg(j,3) = max(y3(end,9));
    end 
        sensitivity(:,i) = (bg(:,3)- bg(:,2));
end

%% dose response

load fitalpha.mat;load Ki;
 drugList = {'Isoproterenol'; 'Basal' ;'Epinephrine'; 'Norepinephrine'; 'Fenoterol' ;'Formoterol';
'Salbutamol' ;'Terbutaline' ;'Broxaterol' ;'Salmeterol'; 'BRL-37344'; 'CGP-12177';
'Alprenolol' ;'Pindolol'; 'SR 59230A'; 'Atenolol' ;'Carvedilol' ;'Metoprolol'; 'Bisoprolol';
'Propranolol' ;'CGP-20712' ;'ICI-118551'};

 Ki = [0.224,NaN,3.97,3.57,13.6,1.71,2.44,31.3,1.31,1.6,37.9,4.70e-03,...
      5.80e-03,2.60e-03,1.64e-02,3.88e-01,1.70e-03,4.70e-02,2.24e-02,...
      1.80e-03,4.50e-03,4.95e-02];
  KLcalc =Ki.*(fitalpha*KR + 1)./(fitalpha*(KR + 1));
   KLcalc = (0.2/KLcalc(1))*KLcalc; 
KG = 2.4131;gamma_L =  0.3762;
p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);

RelTol = 1e-13;MaxStep = 1e-1;
options = odeset('MaxStep',1e3,'NonNegative',[1:29],'RelTol',RelTol);
 y0 = zeros(29,1); y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;
t = [0;1000*60*1000]; p(1:2) = 0;[~,y] = ode15s(@daeODE,t,y0,options,p);y0 = y(end,:);
tspan = [0;2*60*1000]; 
dose = 10.^[-5:.1:3];

beta_blocker = [20 18 17];

for i=1:length(beta_blocker)
    for j = 1:length(dose)
        disp(['Beta Blocker:' drugList{beta_blocker(i)} ' dose: ' num2str(dose(j))]);
        KL =  KLcalc(1);
        alpha_L = fitalpha(1);
        KA = KLcalc(beta_blocker(i)); alpha_A = fitalpha(beta_blocker(i));
        p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A);
        tspan = [0;2*60*1000];
        p(1) = 0; p(2) = 0;  % agonist conc and beta blocker conc
        [t1,y] = ode15s(@daeODE,tspan,y0,options,p);
        bg(j,1) = max(y(end,9));
        yfinal = y(end,:); 
        tspan = [2*60*1000;6*60*1000];
        p(1) = 0.1; p(2) = dose(j);  % agonist conc and beta blocker conc
        [t2,y2] = ode15s(@daeODE,tspan,yfinal,options,p); 
        bg(j,2) = max(y2(end,9));
        yfinal = y2(end,:); 
    end 
        response(:,i) = bg(:,2);
         

end
%% model-experimental comparison

close all;
cd('B:\Robert\Dissertation\Beta Blocker modeling\final_model validation\data\modelKL0p2');
load('ISOPRO_Ca');load('ISOPRO_cAMP');cAMPPRO = cAMP;
load('ISOMET_Ca');load('ISOMET_cAMP');cAMPMET= cAMP;
load('ISOCAR_Ca');load('ISOCAR_cAMP');cAMPCAR = cAMP;
concentration = 10.^[-5:.25:2.5];
for i = 1:length(concentration)
    temp = ISOPRO_Ca{i};temp = (temp - min(temp))./min(temp);temp = findpeaks(temp);temp = temp./temp(1);
    calciumPRO(i) = temp(end);
    temp = ISOMET_Ca{i};temp = (temp - min(temp))./min(temp);temp = findpeaks(temp);temp = temp./temp(1);
    calciumMET(i) = temp(end);
     temp = ISOCAR_Ca{i};temp = (temp - min(temp))./min(temp);temp = findpeaks(temp);temp = temp./temp(1);
    calciumCAR(i) = temp(end);
end
figure;subplot(1,3,1);semilogx(concentration*1e-6,(calciumPRO-min(calciumPRO))./(max(calciumPRO)-min(calciumPRO)));hold all;
subplot(1,3,2);semilogx(concentration*1e-6,(calciumMET-min(calciumMET))./(max(calciumMET)-min(calciumMET)));hold all;
subplot(1,3,3);semilogx(concentration*1e-6,(calciumCAR-min(calciumCAR))./(max(calciumCAR)-min(calciumCAR)));hold all;

cd('B:\Robert\Dissertation\Beta Blocker modeling\final_model validation\data');
load('PRO');PRO(:,2) = (PRO(:,2)-PRO(end,2))./(PRO(2,2)-PRO(end,2));
PRO(:,3) = (PRO(:,3)-PRO(end,2))./(PRO(1,2)-PRO(end,2));PRO(1,1) = concentration(1)*1e-6;
subplot(1,3,1);errorbar(PRO(:,1),PRO(:,2),PRO(:,3),'o','color','g');

load('MET');MET(:,2) = (MET(:,2)-MET(end,2))./(MET(1,2)-MET(end,2));
MET(:,3) = (MET(:,3)-MET(end,2))./(MET(1,2)-MET(end,2));MET(1,1) = concentration(1)*1e-6;
subplot(1,3,2);errorbar(MET(:,1),MET(:,2),MET(:,3),'o','color','g');

load('CAR');CAR(:,2) = (CAR(:,2)-CAR(end,2))./(CAR(1,2)-CAR(end,2));
CAR(:,3) = (CAR(:,3)-CAR(end,2))./(CAR(1,2)-CAR(end,2));CAR(1,1) = concentration(1)*1e-6;
subplot(1,3,3);errorbar(CAR(:,1),CAR(:,2),CAR(:,3),'o','color','g');

%% carvedilol new dose

%run('B:\Robert\Dissertation\Beta Blocker modeling\final_model validation\data\CARsensexp');close all;


 p = daePARAMS(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,scaling_factor);p(1) = 0;
 y0 = zeros(29,1);y0(1) = 0.0132;y0(11) = PKAItot ; y0(17) = PKAIItot;flag = 0;
 [~,y] = ode15s(@daeODE,[0;60*60*1000],y0,[],p,flag);y0Sig = y(end,:);
options = odeset('RelTol',1e-5,'MaxStep',5e-3,'Stats','on');

load y0full;y0 = y0full;
dose = 0.3;drugs = 17; 
KL = KLcalc(1);
alpha_L = fitalpha(1);
KA = KLcalc(drugs); alpha_A = fitalpha(drugs);
p = daePARAMSc(KR,KL,KA,KG,alpha_L,alpha_A,gamma_L,gamma_A,scaling_factor);
    
tspan = [0;2*60];
p(1) = 0; p(2) = 0; flag = 0; % agonist conc and beta blocker conc
[t,y] = ode15s(@daeODEc,tspan,y0,options,p,flag);

bgcAMP{1,1} = y(:,9);
bgCa{1,1} = y(:,47);
time{1,1} = t;
yfinal = y(end,:); 

tspan = [2*60;6*60];
p(1) = 0.1; p(2) = dose;flag = 1;  % agonist conc and beta blocker conc
[t2,y2] = ode15s(@daeODEc,tspan,yfinal,options,p,flag); 
bgcAMP{1,2} = y2(:,9);
bgCa{1,2} = y2(:,47);
time{1,2} = t2;
yfinal = y2(end,:); 
tspan2 = [6*60 10*60]; %changed the time
p(1) = 10; p(2) = dose; flag = 1;    % agonist conc and beta blocker conc
[t3,y3]=ode15s(@daeODEc,tspan2,yfinal,options,p,flag);
bgcAMP{1,3} = y3(:,9);
bgCa{1,3} = y3(:,47);
time{1,3} = t3;


CARmodel = [bgCa{1,1}; bgCa{1,2} ;bgCa{1,3}];CARmodelt = [time{1,1}; time{1,2} ;time{1,3}]./60;

%%
figure;subplot(2,1,1);plot(CARmodelt,CARmodel);box off;xlabel('time (min)');
ylabel('Ca^{2+}')
figure;subplot(2,1,2);
errorbar(1:10,mean(carAa),ste(carAa),'MarkerSize',6,'Marker','o','LineWidth',2); 
ylim([0.49 2.2]);xlim([0.9 10.1]);box off;xlabel('time (min)');ylabel(['Ca^{2+} amp.',sprintf('\n'),' (fold change)']);
set(gcf, 'PaperPositionMode', 'manual');
 set(gcf, 'PaperUnits', 'inches');
 set(gcf, 'PaperPosition', [2.0 4.5 3.8  2.3]);

