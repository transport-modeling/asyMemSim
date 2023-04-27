function [sysInfo] = expSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%specify single simulation
 
%mixture components, compositions, system parameters
  %SBAD1 Data
    sysInfo.memID = "SBAD1";  % polymer spec
%     %9 COMP M1
    sysInfo.mixID = ["TOL","MCH","MNP","DEC","NOC","IOC","TBB","TPB","ICT"];  % mixture spec
    sysInfo.yf = [0.171;0.281;0.0199;0.107;0.221;0.15;0.0217;0.0158;0.013];  % composition spec (must sum to 1)
    PuM = 40; %bar
  %system specs
    sysInfo.n = length(sysInfo.yf);  % number of permeants
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
