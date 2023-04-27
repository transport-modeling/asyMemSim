function [sysInfo] = expSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%specify experiment

%mixture components, compositions, system parameters
  %SBAD1 Data
%     sysInfo.memID = "SBAD1";  % polymer spec
%     sysInfo.mixID_OG = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec **NUMBER of entries should equal params.compID number listed in dataBank.m
%     sysInfo.memID = "PIM1";  % polymer spec
%     sysInfo.mixID_OG = ["TOL","HEP","PXY","OXY","ICT"];
    sysInfo.memID = "MATRIMID";
    sysInfo.mixID_OG = ["TOL","MES","TPB"];
    sysInfo.pressureVec = [10,20,30,40,50,60,70,80,90,100]; %bar -- specify pressures to loop over
    sysInfo.yf = 1;
    
  %system specs
    sysInfo.n = 1;  % number of permeants
    sysInfo.Pd = 1;  % support side pressure (atm)
    sysInfo.pervapMode = 0; % %BETA(UNSTABLE) pervap capability -- yes if == 1 (assume downstream pressure = 0 bar)
    sysInfo.T = 295;  % system temperature (K)
    sysInfo.R = 82.05;  % gas constant (atm cm^3/ mol K)
    sysInfo.diffFit = 0; %see code below for diffusivity fitting
    sysInfo.crossDiffFudge = 0; %if = 0 then no specified cross-diffusivities D_ij will be used
    if sysInfo.crossDiffFudge == 1
        sysInfo.crossDiffSpecs = [1,2,1,3]; % e.g. [i,j] = D_ij^MS
        sysInfo.crossDiffVals = [0.1,0.5]; %um2/s
    end
    sysInfo.lowDiffErrorBar = 0; % if = 1 then diffusivities fit to low error sing comp will be used
%------------------------------------------------------------------------------------------------------------------------------------%
end
