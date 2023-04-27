%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% script:    singCompLoop_asyMemLocal.m                                   %
% Description: Asymmetric membrane local flux model/method solver for:    %
%                -single component permeation predictions loop            %
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
%pure comp predictions loop
    pressureVec = sysInfo.pressureVec; %bar
    modSpec = sysInfo.modSpec;
    [allMPSCF,failedRuns] = asyMemSingCompEval(sysInfo,pressureVec,modSpec); %MPSCF = model, pressure, single component fluxes
%------------------------------------------------------------------------------------------------------------------------------------%

    
    