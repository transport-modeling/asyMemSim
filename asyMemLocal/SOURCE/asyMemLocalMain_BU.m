%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script:    asyMemLocalMain.m                                            %
% Description: Asymmetric membrane local flux model/method solver for:    %
%                -single mixture permeation simulation                    %
%                -single component diffsivity fitting from single exp     %
%                -single component permeation predictions loop            %
%                -mixture permeation prediciotns loop                     %
%                (See tutorial document for specific info regarding use)  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------------------------------------------------------------%
%specify single simulation
 
%mixture components, compositions, system parameters
  %SBAD1 Data
    sysInfo.memID = "SBAD1";  % polymer spec
    %9 COMP M1
    sysInfo.mixID = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
    sysInfo.yf = [0.171;0.281;0.0199;0.107;0.221;0.15;0.0217;0.0158;0.013];  % composition spec (must sum to 1)
    PuM = 40; %bar
    %5 COMP M1
%     sysInfo.mixID = ["TOL","MNP","NOC","IOC","ICT"];  % mixture spec
%     sysInfo.yf = [0.24;0.0239;0.453;0.266;0.017];  % composition spec (must sum to 1)   
%     PuM = 30; %bar
%     sysInfo.crossDiffAnal = 1;
%     sysInfo.crossDiffSpecs = [1,3,1,4]; % e.g. [i,j] = D_ij^MS
%     sysInfo.crossDiffVals = [0.02,0.2]; %um2/s
    %3 COMP M1
    sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
    sysInfo.yf = [0.282;0.385;0.333];  % composition spec (must sum to 1)
%     PuM = 30; %bar
    PuM = 60; %bar
%     sysInfo.crossDiffFudge = 1;
%     sysInfo.crossDiffSpecs = [1,2,1,3]; % e.g. [i,j] = D_ij^MS
%     sysInfo.crossDiffVals = [0.1,0.5]; %um2/s
    %3 COMP M2
%     sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
%     sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
% %     sysInfo.yf = [0.531;0.363;0.104];  % composition spec (must sum to 1)
%     sysInfo.yf = [0.748;0.224;0.028];
%     sysInfo.yf = [0.763;0.197;0.040];
%      sysInfo.yf = [0.766;0.196;0.038];
%     Pum = 1;
%       PuM = 15; %bar
%   PuM = 30; %bar
    %3 COMP M3
%     sysInfo.mixID = ["TOL","NOC","ICT"];  % mixture spec
%     sysInfo.yf = [0.04;0.95;0.01];  % composition spec (must sum to 1)
% %     sysInfo.yf = [0.976;0.01514;0.0089];  % composition spec (must sum to 1)
%     PuM = 30; %bar
    %2 COMP M1
%     sysInfo.mixID = ["TOL","TPB"];  % mixture spec
%     sysInfo.yf = [0.9904;0.0096];  % composition spec (must sum to 1)
%     PuM = 30; %bar

  %PIM1 Data
%     sysInfo.memID = "PIM1";  % polymer spec
    %5 COMP M1
%     sysInfo.mixID = ["TOL","HEP","PXY","OXY","ICT"]; % mixture spec
%     sysInfo.yf = [0.257;0.216;0.205;0.264;0.058];  % composition spec (must sum to 1)
%     PuM = 30; %bar 
%     %3 COMP M1
%     sysInfo.mixID = ["TOL","PXY","ICT"]; % mixture spec
%     sysInfo.yf = [0.406;0.401;0.193]; % composition spec (must sum to 1)
%     PuM = 1; %bar
%     %3 COMP M2
%     sysInfo.mixID = ["TOL","HEP","PXY"]; % mixture spec
%     sysInfo.yf = [0.345;0.362;0.293]; % composition spec (must sum to 1)
%     PuM = 1; %bar
      %2 COMP M1
%     sysInfo.mixID = ["HEP","OXY"]; % mixture spec
%     sysInfo.yf = [0.2656;0.7341]; % composition spec (must sum to 1)
%     sysInfo.yf = [0.4866;0.5134];
%     sysInfo.yf = [0.7682;0.2318];
%     PuM = 1; %bar


  %single component selection   
%     sysInfo.mixID = ["HEP"];
%     PuM = 60; %bar
%     sysInfo.yf = [1];  % composition spec (must sum to 1)

  %system specs
    sysInfo.n = length(sysInfo.yf);  % number of permeants
    sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
    sysInfo.Pd = 1;  % support side pressure (atm)
    sysInfo.T = 295;  % system temperature (K)
    sysInfo.R = 82.05;  % gas constant (atm cm^3/ mol K)
    sysInfo.diffFit = 0; %see code below for diffusivity fitting
    sysInfo.crossDiffFudge = 0; %if = 0 then no specified cross-diffusivities D_ij will be used
    sysInfo.lowDiffErrorBar = 0; % if = 1 then diffusivities fit to low error sing comp will be used
        
  %sorption membrane phase model
    %sysInfo.memPhaseModel = "F-H";
    %sysInfo.memPhaseModel = "DSM";
    sysInfo.memPhaseModel = "FH-LM";
    
  %specify thermodynamic assumptions
    sysInfo.noThermoCoupling = 0;
    %sysInfo.noThermoCoupling = 1;
    sysInfo.memPhaseModel_FicksOG = 0;
    %sysInfo.memPhaseModel_FicksOG = 1;
    
  %specify diffusional relationships    
    sysInfo.diffModel = "NoCoupling";  
    %sysInfo.diffModel = "Vignes";
    %sysInfo.diffModel = "Darken";
    %sysInfo.diffModel = "Fudge";  
    %sysInfo.swlDiffModel = "none";
    %sysInfo.swlDiffModel = "FFV-OG";
    %sysInfo.swlDiffModel = "FFV-UAV";
    sysInfo.swlDiffModel = "Avg-Diff";
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%cross-diffusion analysis
%     sysInfo.crossDiffFudge = 1;
%     sysInfo.crossDiffSpecs = [1,2,1,3]; % e.g. [i,j] = D_ij^MS
%     sysInfo.crossDiffVals = [0.1,0.5]; %um2/s
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%numerical method sepcs
  %shooting algorithm (SA) [Recommended]
    sysInfo.numMethod = "ShootAlg"; 
  %full-discretization (FUD)
%     sysInfo.numMethod = "FullDis";
%     sysInfo.numNodes = 20; 
  %multiple shooting (MSA)
%     sysInfo.numMethod = "MultShootAlg";
%     sysInfo.numShootPoints = 1;  
%     sysInfo.casADi = 1; % 0 = fsolve
  %solver specs
    sysInfo.solverSpec = "levenberg-marquardt"; %default, use for FUD
%     sysInfo.solverSpec = "trust-region-dogleg"; %use if other solver fails
    sysInfo.iterDetail = 1; %0 = all main solver output will be suppressed 
  %initial guess specs
    sysInfo.initGuessApprox = 1; %0 = use feed mol fraction vector and small total flux value (use if IVP approx fails)
%------------------------------------------------------------------------------------------------------------------------------------%    

%------------------------------------------------------------------------------------------------------------------------------------%    
%run single simulation
%       [localCompFlux,partialFlux] = asyMemLocalSolve(sysInfo)
%------------------------------------------------------------------------------------------------------------------------------------%    

%------------------------------------------------------------------------------------------------------------------------------------%
%diffusivity fit from sing comp permeation
    sysInfo.initGuessApprox = 0;
    sysInfo.iterDetail = 0; %0 = all main solver output will be suppressed
    sysInfo.memPhaseModel_FicksOG = 1;
    for i = 1:2
        if i == 1
%           SABD1
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID_OG = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
            singFluxVec = [8.76;0.74;0.054;0.063;2.23;0.302;0.23;0.0054;0.0266]; %LMH 20 bar flux (30 IOCT)
            singFluxVec_Error = [0.74;0.36;0.011;0.014;0.5;0.120;0.07;0;0.0024]; %LMH 20 bar flux error (30 IOCT)
%             singFluxVec = [27.53;2.37;0.121;0.186;7.5;0.627;0.66;0.0233;0.0762]; %LMH 60 bar flux (30 IOCT)
%             singFluxVec_Error = [3.09;0.93;0.008;0.015;1.71;0.3;0.16;0;0.0037]; %LMH 60 bar flux error (30 IOCT)
            singFluxVec_LOW = [singFluxVec-singFluxVec_Error]; %LMH
%           run fitting function
            modSpec = [3]; %FH - 1 ; DSM - 2; FH-LM - 3
            [allModSingCompDiff_SBAD,allModSingCompDiffLOW_SBAD] = ...
                asyMemDiffFitSolve(sysInfo,singFluxVec,singFluxVec_Error,modSpec);
        elseif i == 2
%           PIM1
            sysInfo.memID = "PIM1";  % polymer spec
            sysInfo.mixID_OG = ["TOL","HEP","PXY","OXY","ICT"];
            singFluxVec = [10.04;27.7;7.88;4.74;0.93]; %LMH %20baru flux
            singFluxVec_Error = [2.58;5.20;1.98;0.87;0.424]; %LMH %20bar flux error
%           run fitting function
            modSpec = [3]; %FH - 1 ; DSM - 2; FH-LM - 3
            [allModSingCompDiff_PIM,allModSingCompDiffLOW_PIM] = ...
                asyMemDiffFitSolve(sysInfo,singFluxVec,singFluxVec_Error,modSpec);
        end
    end
%------------------------------------------------------------------------------------------------------------------------------------%
 
%------------------------------------------------------------------------------------------------------------------------------------%
%mixture simulations and single comp predictions
  %pure comp predictions
%     sysInfo.initGuessApprox = 0;
%     sysInfo.iterDetail = 0; %0 = all main solver output will be suppressed 
%     sysInfo.n = 1;
%     sysInfo.yf = [1];
%     pressureVec = [20,30,40,50,60]; %bar
%     modSpec = [1 2 3];  %FH - 1 ; DSM - 2; FH-LM - 3
%     for i = 1:2
%         if i == 1
%           %SABD1
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID_OG = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
%           %MS Eval
%             sysInfo.memPhaseModel_FicksOG = 0;
%             [allMPSCF_SBAD_MS] = asyMemSingCompEval(sysInfo,pressureVec,modSpec);
%           %Ficks Eval
%             sysInfo.memPhaseModel_FicksOG = 1;
%             [allMPSCF_SBAD_FICKS] = asyMemSingCompEval(sysInfo,pressureVec,modSpec);
%         elseif i == 2
%           %PIM1
%             sysInfo.memID = "PIM1";  % polymer spec
%             sysInfo.mixID_OG = ["TOL","HEP","PXY","OXY","ICT"];
%           %MS Eval
%             sysInfo.memPhaseModel_FicksOG = 0;
%             [allMPSCF_PIM_MS] = asyMemSingCompEval(sysInfo,pressureVec,modSpec);
%           %Ficks Eval
%             sysInfo.memPhaseModel_FicksOG = 1;
%             [allMPSCF_PIM_FICKS] = asyMemSingCompEval(sysInfo,pressureVec,modSpec);
%         end
%     end

  %mixture prediction loop
%     sysInfo.initGuessApprox = 0;
%     sysInfo.iterDetail = 0; %0 = all main solver output will be suppressed 
%     allMixDiffSorpModSim = [];
%     for i = 1:12 %mixture 1,2,3,etc.
%         if i == 1
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
%             sysInfo.yf = [0.171;0.281;0.0199;0.107;0.221;0.15;0.0217;0.0158;0.013];  % composition spec (must sum to 1)
%             PuM = 40; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 2
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
%             sysInfo.yf = [0.284;0.388;0.328];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;       
%         elseif i == 3
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
%             sysInfo.yf = [0.282;0.385;0.333];  % composition spec (must sum to 1)
%             PuM = 60; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 4
%             sysInfo.memID = "PIM1";  % polymer spec
%             sysInfo.mixID = ["TOL","HEP","PXY","OXY","ICT"]; % mixture spec
%             sysInfo.yf = [0.257;0.216;0.205;0.264;0.058];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 5
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","MNP","NOC","IOC","ICT"];  % mixture spec
%             sysInfo.yf = [0.24;0.0239;0.453;0.266;0.017];  % composition spec (must sum to 1)   
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 6
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
%             sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
%             PuM = 15; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 7
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
%             sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 8
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
%             sysInfo.yf = [0.532;0.363;0.105];  % composition spec (must sum to 1)
%             PuM = 15; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 9
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
%             sysInfo.yf = [0.533;0.363;0.104];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 10
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","ICT"];  % mixture spec
%             sysInfo.yf = [0.04;0.95;0.01];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         elseif i == 11
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","NOC","ICT"];  % mixture spec
%             sysInfo.yf = [0.976;0.01514;0.0089];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         else
%             sysInfo.memID = "SBAD1";  % polymer spec
%             sysInfo.mixID = ["TOL","TPB"];  % mixture spec
%             sysInfo.yf = [0.9904;0.0096];  % composition spec (must sum to 1)
%             PuM = 30; %bar
%             sysInfo.crossDiffAnal = 0;
%         end
%         sysInfo.n = length(sysInfo.yf);  % number of permeants
%         sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
%         allTA = [];
%         for j = 1:3 %full MS, Ficks, noThermoCoup
%             if j == 1
%                 sysInfo.noThermoCoupling = 0;
%                 sysInfo.memPhaseModel_FicksOG = 0;
%             elseif j == 2
%                 sysInfo.memPhaseModel_FicksOG = 1;
%             else
%                 sysInfo.noThermoCoupling = 1;
%                 sysInfo.memPhaseModel_FicksOG = 0;
%             end
%             allSM = [];
%             for l = 1:3 %FH, DSM, FH-LM
%                 if l == 1
%                     sysInfo.memPhaseModel = "F-H";
%                 elseif l == 2
%                     sysInfo.memPhaseModel = "DSM";
%                 else
%                     sysInfo.memPhaseModel = "FH-LM";
%                 end
%                 allCrD = [];
%                 for k = 1:2 %Vignes or noDiffCoup
%                     if j == 2 && k == 1
%                         continue
%                     else
%                         if k == 1
%                             sysInfo.diffModel = "Vignes";
%                         else
%                             sysInfo.diffModel = "Fudge";  
%                         end
%                         allCD = [];
%                         for m = 2 %FFV-UAV, D_avg, or no D_im model
%                             if m == 1
%                                 sysInfo.swlDiffModel = "FFV-UAV";
%                             elseif m == 2
%                                 sysInfo.swlDiffModel = "Avg-Diff";
%                             else
%                                 sysInfo.swlDiffModel = "none";
%                             end 
%                             [localCompFlux_m,partialFlux_m] = asyMemLocalSolve(sysInfo);
%                             allCD = [allCD,localCompFlux_m,[partialFlux_m;0]];
%                         end
%                         allCrD = [allCrD,allCD];
%                     end    
%                 end
%                allSM = [allSM,allCrD]; 
%             end
%             allTA = [allTA,allSM];
%         end
%         allMixDiffSorpModSim = [allMixDiffSorpModSim;allTA];
%     end
%------------------------------------------------------------------------------------------------------------------------------------%
    
    