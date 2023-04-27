function [sysInfo] = expSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%specify exp
 
%mixture components, compositions, system parameters
  %SBAD1 Data
%     sysInfo.memID = "SBAD1";  % polymer spec
%     sysInfo.mixID_OG = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec **NUMBER of entries should equal params.compID number listed in dataBank.m
%     sysInfo.singFluxVec = [8.76;0.74;0.054;0.063;2.23;0.302;0.23;0.0054;0.0266]; %LMH 20 bar flux (30 IOCT)
%     sysInfo.singFluxVec_Error = [0.74;0.36;0.011;0.014;0.5;0.120;0.07;0;0.0024]; %LMH 20 bar flux error (30 IOCT)
%     sysInfo.memID = "PIM1";  % polymer spec
%     sysInfo.mixID_OG = ["TOL","HEP","PXY","OXY","ICT"];
%     sysInfo.singFluxVec = [10.04;27.7;7.88;4.74;0.93]; %LMH %20baru flux
%     sysInfo.singFluxVec_Error = [2.58;5.20;1.98;0.87;0.424]; %LMH %20bar flux error
%     PuM = 20; %feed-side pressure (bar)
    sysInfo.memID = "MATRIMID";  % polymer spec
    sysInfo.mixID_OG = ["TOL","MES","TPB"];
    sysInfo.singFluxVec = [0.798;0.0666;0.00734]; %LMH %20baru flux
    sysInfo.singFluxVec_Error = [0;0;0]; %LMH %20bar flux error
    PuM = 20; %feed-side pressure (bar)
  %system specs
    sysInfo.n = 1;  % number of permeants
    sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
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
