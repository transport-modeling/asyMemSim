%TEST SCRIPT LOOP MODEL FOR SIMULATED DATA GENERATION  

%------------------------------------------------------------------------------------------------------------------------------------%
%specify single simulation
 
%mixture components, compositions, system parameters
  %SBAD1 Data
    sysInfo.memID = 'SBAD1';  % polymer spec
    sysInfo.lmem = .3;  % thickness of active membrane layer (um) 
    %9 COMP
    sysInfo.mixID = ['TOL','MCH','MNP','DEC','NOC','IOC','TBB','TPB','ICT'];  % mixture spec
    sysInfo.yf = [0.171;0.281;0.0199;0.107;0.221;0.15;0.0217;0.0158;0.013];  % composition spec (must sum to 1)
    PuM = 40; %bar
    %5 COMP M1
%     sysInfo.mixID = ['TOL','MNP','NOC','IOC','ICT'];  % mixture spec
%     sysInfo.yf = [0.24;0.0239;0.453;0.266;0.017];  % composition spec (must sum to 1)   
%     PuM = 30; %bar
    %3 COMP M1
%     sysInfo.mixID = ['TOL','IOC','ICT'];  % mixture spec
%     sysInfo.yf = [0.282;0.385;0.333];  % composition spec (must sum to 1)
%     PuM = 30; %bar
%     PuM = 60; %bar
    %3 COMP M2
%     sysInfo.mixID = ['TOL','NOC','MNP'];  % mixture spec
%     sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
%     sysInfo.yf = [0.531;0.363;0.104];  % composition spec (must sum to 1)
%     sysInfo.yf = [0.748;0.224;0.028];
%     Pum = 1;
%     PuM = 15; %bar
%     PuM = 30; %bar

  %PIM1 Data
%     sysInfo.memID = 'PIM1';  % polymer spec
%     sysInfo.lmem = 1.5;  % thickness of active membrane layer (um) 
    %5 COMP M1
%     sysInfo.mixID = ['TOL','HEP','PXY','OXY','ICT'];
%     sysInfo.yf = [0.257;0.216;0.205;0.264;0.058];  % composition spec (must sum to 1)
%     PuM = 30; %bar 
    %3 COMP M1
%     sysInfo.mixID = ['TOL','PXY','ICT'];
%     sysInfo.yf = [0.406;0.401;0.193];
%     Pum = 1;

  %single component selection   
%     sysInfo.mixID = ['NOC'];
%     PuM = 20;
%     sysInfo.yf = [1];  % composition spec (must sum to 1)

  %system specs
    sysInfo.n = length(sysInfo.yf);  % number of permeants
    sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
    sysInfo.Pd = 1;  % support side pressure (atm)
    sysInfo.T = 298;  % system temperature (K)
    sysInfo.R = 82.05;  % gas constant (atm cm^3/ mol K)
    sysInfo.diffFit = 0; %see code below for diffusivity fitting
    sysInfo.crossDiffAnal = 0;
        
  %sportion membrane phase model
    %sysInfo.memPhaseModel = 'F-H';
    %sysInfo.memPhaseModel = 'DSM';
    sysInfo.memPhaseModel = 'FH-LM';
    
  %specify thermodynamic assumptions
    sysInfo.noThermoCoupling = 0;
    %sysInfo.noThermoCoupling = 1;
    %sysInfo.memPhaseModel_FicksOG = 0;
    sysInfo.memPhaseModel_FicksOG = 1;
    
  %specify diffusional relationships    
    sysInfo.diffModel = 'NoCoupling';  
    %sysInfo.diffModel = 'Vignes';
    %sysInfo.diffModel = 'Darken';
    %sysInfo.diffModel = 'Fudge';  
    sysInfo.swlDiffModel = 'none';
    %sysInfo.swlDiffModel = 'FFV-OG';
    %sysInfo.swlDiffModel = 'FFV-UAV';
    %sysInfo.swlDiffModel = 'Avg-Diff';
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%cross-diffusion analysis
%     sysInfo.crossDiffAnal = 1;
%     sysInfo.crossDiffSpecs = [2,5,3,5,4,5]; % e.g. [i,j] = D_ij^MS
%     sysInfo.crossDiffVals = [0.9,0.09,0.001]; %um2/s
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%    
%run single simulation
      [localCompFlux,partialFlux] = asyMemLocalSolve(sysInfo)
%------------------------------------------------------------------------------------------------------------------------------------%    

%------------------------------------------------------------------------------------------------------------------------------------%
%isothem generator
%     compVec = [];
%     for i = 1:20
%         CompVec = [0.01;0.05;0.1;0.15;0.2;0.25;0.3;0.35;0.4;0.45;0.5;0.55;0.6;0.65;0.7;0.75;0.8;0.85;0.9;0.95;0.99];
%         sysInfo.yf = [CompVec(i);1-CompVec(i)];
%         [~,~,phiFeed] = asyMemLocalSolve(sysInfo.);
%         compVec = [compVec,phiFeed];
%     end
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%diffusivity fit from sing comp permeation
%   %SABD1
% %     sysInfo.memID = 'SBAD1';  % polymer spec
% %     sysInfo.lmem = .3;  % thickness of active membrane layer (um)
% %     sysInfo.mixID_OG = ['TOL','MCH','MNP','DEC','NOC','IOC','TBB','TPB','ICT'];  % mixture spec
% %     singFluxVec = [8.76;0.74;0.054;0.063;2.23;0.302;0.23;0.0054;0.0266]; %LMH 20 bar flux (30 IOCT)
% %     singFluxVec_Error = [0.74;0.36;0.011;0.014;0.5;0.120;0.07;0;0.0024]; %LMH 20 bar flux error (30 IOCT)
% % %     singFluxVec = [27.53;2.37;0.121;0.186;7.5;0.627;0.66;0.0233;0.0762]; %LMH 60 bar flux (30 IOCT)
% % %     singFluxVec_Error = [3.09;0.93;0.008;0.015;1.71;0.3;0.16;0;0.0037]; %LMH 60 bar flux error (30 IOCT)
% %     singFluxVec_LOW = [singFluxVec-singFluxVec_Error]; %LMH
%   %PIM1
%     sysInfo.memID = 'PIM1';  % polymer spec
%     sysInfo.lmem = 1.5;  % thickness of active membrane layer (um) 
%     sysInfo.mixID_OG = ['TOL','HEP','PXY','OXY','ICT'];
%     singFluxVec = [10.04;27.7;7.88;4.74;0.93]; %LMH %20baru flux
%     singFluxVec_Error = [2.58;5.20;1.98;0.87;0.424]; %LMH %20bar flux error
%   %run fitting function
%     modSpec = [1 2 3]; %FH - 1 ; DSM - 2; FH-LM - 3
%     [allModSingCompDiff,allModSingCompDiffLOW] = ...
%     asyMemDiffFitSolve(sysInfo,singFluxVec,singFluxVec_Error,modSpec);
%------------------------------------------------------------------------------------------------------------------------------------%
 
%------------------------------------------------------------------------------------------------------------------------------------%
%mixture simulations and single comp predictions
  %pure comp predictions
%     sysInfo.n = 1;
%     sysInfo.yf = [1];
%     %SABD1
% %     sysInfo.memID = 'SBAD1';  % polymer spec
% %     sysInfo.mixID_OG = ['TOL','MCH','MNP','DEC','NOC','IOC','TBB','TPB','ICT'];  % mixture spec
%     %PIM1
%     sysInfo.memID = 'PIM1';  % polymer spec
%     sysInfo.mixID_OG = ['TOL','HEP','PXY','OXY','ICT'];
%   
%     pressureVec = [20,30,40,50,60];
%     modSpec = [1 2 3];
%     [allMPSCF] = asyMemSingCompEval(sysInfo,pressureVec,modSpec);

%------------------------------------------------------------------------------------------------------------------------------------%
    
    