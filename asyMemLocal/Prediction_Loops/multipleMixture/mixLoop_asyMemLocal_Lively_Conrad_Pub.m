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
    for i = [1 2] %mixture 1,2,3,etc.
        if i == 1
            sysInfo.memID = "MATRIMID";  % polymer spec
            sysInfo.mixID = ["TOL","MES"];  % mixture spec
            sysInfo.yf = [0.98;0.02];  % composition spec (must sum to 1)
            sysInfo.crossDiffFudge = 0;
        elseif i == 2
            sysInfo.memID = "MATRIMID";  % polymer spec
            sysInfo.mixID = ["TOL","TPB"];  % mixture spec
            sysInfo.yf = [0.98;0.02];  % composition spec (must sum to 1)
            sysInfo.crossDiffFudge = 0; 
        elseif i == 3
            sysInfo.memID = "MATRIMID";  % polymer spec
            sysInfo.mixID = ["TOL"];  % mixture spec
            sysInfo.yf = [1];  % composition spec (must sum to 1)
            sysInfo.crossDiffFudge = 0; 
        elseif i == 4
            sysInfo.memID = "MATRIMID";  % polymer spec
            sysInfo.mixID = ["MES"];  % mixture spec
            sysInfo.yf = [1];  % composition spec (must sum to 1)
            sysInfo.crossDiffFudge = 0; 
        else 
            sysInfo.memID = "MATRIMID";  % polymer spec
            sysInfo.mixID = ["TPB"];  % mixture spec
            sysInfo.yf = [1];  % composition spec (must sum to 1)
            sysInfo.crossDiffFudge = 0; 
        end
        sysInfo.n = length(sysInfo.yf);  % number of permeants
        allTA = [];
        for j = 1:10 %full MS, Ficks, noThermoCoup
            PuM = [10 20 30 40 50 60 70 80 90 100]; 
            sysInfo.Pu = PuM(j)*0.9869;  % feed side pressure [bar]*0.9832 = [atm]            
            allSM = [];
            sysInfo.noThermoCoupling = 0;
            sysInfo.memPhaseModel_FicksOG = 0;
            sysInfo.noGammaFugacityODEs = 0;
            for l = 1 %FH, DSM, FH-LM
                if l == 1
                    sysInfo.memPhaseModel = "F-H";
                elseif l == 2
                    sysInfo.memPhaseModel = "DSM";
                else
                    sysInfo.memPhaseModel = "FH-LM";
                end
                allCrD = [];
                for k = 1 %Vignes or noDiffCoup
                    if j == 100 && k == 1
                        continue
                    else
                        if k == 1
                            sysInfo.diffModel = "Vignes";
                        else
                            sysInfo.diffModel = "NoCoupling";  
                        end
                        allCD = [];
                        for m = 3 %FFV-UAV, D_avg, or no D_im model
                            if m == 1
                                sysInfo.swlDiffModel = "FFV";
                            elseif m == 2
                                sysInfo.swlDiffModel = "Avg-Diff";
                            else
                                sysInfo.swlDiffModel = "none";
                            end 
                            [localCompFlux_m,partialFlux_m,exitFlag,params] = asyMemLocalSolve(sysInfo);
                            purFug_Pd = params.psat.*exp(params.Vs(1:params.n).*(params.Pd-params.psat*0.00131)/(params.R*params.T));
                            mixFug_Pu = params.yf.*params.psat.*exp(params.Vs(1:params.n).*(params.Pu-params.psat*0.00131)/(params.R*params.T));
                            mixFug_Pd = localCompFlux_m(1:2).*(exp(-params.Vs(1:params.n).*(params.Pu-params.Pd)/(params.R*params.T))).*purFug_Pd;
                            phiFH_Pu = y2phiPhaseEq_FH(params);
                            options = optimoptions(@fsolve,'Display','off','FunctionTolerance',1E-6,'MaxFunctionEvaluations',4500,...
                                'MaxIterations',500,'Algorithm',params.solverSpec);
                            phiFH_Pd = fsolve(@(phis)[localCompFlux_m(1:params.n)-phis2yPhaseEq_FH_RHS([phis;1-sum(phis)],params)],[0.5;0.001],options);
                            phiFH_Pd(end+1) = 1-sum(phiFH_Pd);   
                            
                            allCD = [allCD,localCompFlux_m(end),phiFH_Pu(1),phiFH_Pd(1),phiFH_Pu(2),phiFH_Pd(2),mixFug_Pu(1),mixFug_Pd(1),mixFug_Pu(2),mixFug_Pd(2),localCompFlux_m(1)];
%                             allCD = [allCD,localCompFlux_m(end),phiFH_Pu(1),phiFH_Pd(1),mixFug_Pu(1),mixFug_Pd(1)];
                            if exitFlag ~= 1 && exitFlag ~= 2 && exitFlag ~= 3 && exitFlag ~= 4
                                failedRuns = [failedRuns,[exitFlag;i;j;l;k;m]];
                            end                         
                        end
                        allCrD = [allCrD,allCD,];
                    end    
                end
               allSM = [allSM,allCrD]; 
            end
            allTA = [allTA;allSM];
        end
        allMixDiffSorpModSim = [allMixDiffSorpModSim;allTA];
    end
    
        
        
%------------------------------------------------------------------------------------------------------------------------------------%
    
    