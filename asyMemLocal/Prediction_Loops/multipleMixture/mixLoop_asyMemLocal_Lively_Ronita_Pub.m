%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script:    mixLoop_asyMemLocal.m                                        %
% Description: Asymmetric membrane local flux model/method solver for:    %
%                -mixture permeation prediciotns loop                     %
%                (See tutorial document for specific info regarding use)  %
%                                                                         %
% Copyright© 2021 Georgia Tech Research Corporation                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%------------------------------------------------------------------------------------------------------------------------------------%
%add /SOURCE/ folder and /casADi-v3.5.5/ to MATLAB session search path
    soFolder = fullfile(pwd,'..','..','SOURCE');
    mainFolder = fullfile(pwd,'..','..');
    addpath(fullfile(soFolder,'casadi-v3.5.5'));
    addpath(genpath(soFolder));
    addpath(genpath(mainFolder));
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%specify simulation input file names and load data
    sysInfoEXP = expSpec_EXAMPLE(); %change function name to match input file name 
    sysInfoSMS = solverModelSpec_EXAMPLE(); %change function name to match input file name
    sysInfo = cell2struct([struct2cell(sysInfoEXP);struct2cell(sysInfoSMS)],[fieldnames(sysInfoEXP);fieldnames(sysInfoSMS)]);
    sysInfo.dataBankName = "dataBank_EXAMPLE";
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------%
%mixture prediction loop
    allMixDiffSorpModSim = [];
    failedRuns = [];
    for i = [1 2 3 5] %mixture 1,2,3,etc.
        if i == 1
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
            sysInfo.yf = [0.171;0.281;0.0199;0.107;0.221;0.15;0.0217;0.0158;0.013];  % composition spec (must sum to 1)
            PuM = 40; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 2
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
            sysInfo.yf = [0.284;0.388;0.328];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;       
        elseif i == 3
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
            sysInfo.yf = [0.282;0.385;0.333];  % composition spec (must sum to 1)
            PuM = 60; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 4
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","IOC","ICT"];  % mixture spec
            sysInfo.yf = [0.282;0.385;0.333];  % composition spec (must sum to 1)
            PuM = 60; %bar
            sysInfo.crossDiffFudge = 1;
            sysInfo.crossDiffSpecs = [1,2,1,3]; % e.g. [i,j] = D_ij^MS
            sysInfo.crossDiffVals = [0.1,0.5]; %um2/s
        elseif i == 5
            sysInfo.memID = "PIM1";  % polymer spec
            sysInfo.mixID = ["TOL","HEP","PXY","OXY","ICT"]; % mixture spec
            sysInfo.yf = [0.257;0.216;0.205;0.264;0.058];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 6
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","MNP","NOC","IOC","ICT"];  % mixture spec
            sysInfo.yf = [0.24;0.0239;0.453;0.266;0.017];  % composition spec (must sum to 1)   
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 7
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
            sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
            PuM = 15; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 8
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
            sysInfo.yf = [0.793;0.167;0.039];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 9
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
            sysInfo.yf = [0.532;0.363;0.105];  % composition spec (must sum to 1)
            PuM = 15; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 10
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","MNP"];  % mixture spec
            sysInfo.yf = [0.533;0.363;0.104];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 11
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","ICT"];  % mixture spec
            sysInfo.yf = [0.04;0.95;0.01];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        elseif i == 12
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","NOC","ICT"];  % mixture spec
            sysInfo.yf = [0.976;0.01514;0.0089];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        else
            sysInfo.memID = "SBAD1";  % polymer spec
            sysInfo.mixID = ["TOL","TPB"];  % mixture spec
            sysInfo.yf = [0.9904;0.0096];  % composition spec (must sum to 1)
            PuM = 30; %bar
            sysInfo.crossDiffFudge = 0;
        end
        sysInfo.n = length(sysInfo.yf);  % number of permeants
        sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
        allTA = [];
        for j = 1:3 %full MS, Ficks, noThermoCoup
            if j == 1
                sysInfo.noThermoCoupling = 0;
                sysInfo.memPhaseModel_FicksOG = 0;
                sysInfo.noGammaFugacityODEs = 1; % ONLY USE FOR FULL MS (FICK == 0 and noThermo == 0)! option keep MS system interms of fugacity gradients and do not make chain rule trick to get gamma matrix
            elseif j == 2
                sysInfo.memPhaseModel_FicksOG = 1;
                sysInfo.noGammaFugacityODEs = 0;
            else
                sysInfo.noThermoCoupling = 1;
                sysInfo.memPhaseModel_FicksOG = 0;
                sysInfo.noGammaFugacityODEs = 0;
            end
            allSM = [];
            for l = 1:3 %FH, DSM, FH-LM
                if l == 1
                    sysInfo.memPhaseModel = "F-H";
                elseif l == 2
                    sysInfo.memPhaseModel = "DSM";
                else
                    sysInfo.memPhaseModel = "FH-LM";
                end
                allCrD = [];
                for k = 1:2 %Vignes or noDiffCoup
                    if j == 2 && k == 1
                        continue
                    else
                        if k == 1
                            sysInfo.diffModel = "Vignes";
                        else
                            sysInfo.diffModel = "NoCoupling";  
                        end
                        allCD = [];
                        for m = 2 %FFV-UAV, D_avg, or no D_im model
                            if m == 1
                                sysInfo.swlDiffModel = "FFV";
                            elseif m == 2
                                sysInfo.swlDiffModel = "Avg-Diff";
                            else
                                sysInfo.swlDiffModel = "none";
                            end 
                            [localCompFlux_m,partialFlux_m,exitFlag] = asyMemLocalSolve(sysInfo);
                            allCD = [allCD,localCompFlux_m,[partialFlux_m;0]];
                            if exitFlag ~= 1 && exitFlag ~= 2 && exitFlag ~= 3 && exitFlag ~= 4
                                failedRuns = [failedRuns,[exitFlag;i;j;l;k;m]];
                            end
                        end
                        allCrD = [allCrD,allCD];
                    end    
                end
               allSM = [allSM,allCrD]; 
            end
            allTA = [allTA,allSM];
        end
        allMixDiffSorpModSim = [allMixDiffSorpModSim;allTA];
    end
%------------------------------------------------------------------------------------------------------------------------------------%
    
    