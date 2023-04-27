function [sysInfo] = expSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%specify single simulation
 
%mixture components, compositions, system parameters
  %SBAD1 Data
    sysInfo.memID = "Ducky-10-12comp";  % polymer spec
    %12 COMP M1
    sysInfo.mixID = ["IOC","TOL","BCH","TET","MNP","ICT","TPB","PRI","DDB","PYR","DCS","SQU" ];  % mixture spec
    
    %Albatross-5-12comp feed composition below
    %sysInfo.yf = [0.0193;0.8924;0.0064;0.0092;0.0077;0.0117;0.0107;0.0071;0.0084;0.0074;0.0098;0.01];  % composition spec (must sum to 1)
    
    %Ducky-10-12comp feed composition below
    %sysInfo.yf = [0.0126767015392605;0.874955024936112;0.0107799133222996;0.0112794710793953;0.0109978527000903;0.0111994143059899;0.011689929770511;0.0111221962999486;0.0111741546080492;0.0109983150008092;0.0114540441780076;0.0116729822595269];  % composition spec (must sum to 1)
    
    %Ducky-9-12comp feed composition below
    sysInfo.yf = [0.0103276162247238;0.883434986215907;0.0102433116826168;0.0104430839738515;0.0101305616107897;0.0109774377992201;0.011252857686684;0.0107043883954339;0.0104608400849019;0.00927498940308899;0.011060428809073;0.0116894981137093];  % composition spec (must sum to 1)
    
    %Mallard-1-12comp feed composition below
    %sysInfo.yf = [0.01370891122;0.872133954;0.01169827194;0.01152161375;0.0115658596;0.01264385258;0.01316118212;0.00509529692;0.01210352992;0.01168821572;0.01196063483;0.0127186774];  % composition spec (must sum to 1)
    
    PuM = 30; %feed-side pressure (bar)
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
