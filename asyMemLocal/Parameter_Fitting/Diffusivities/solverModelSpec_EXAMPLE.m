function [sysInfo] = solverModelSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%sorption/diffusion model and thermodynamic assumpstions spec         
  %sorption membrane phase model
    sysInfo.modSpec = [1]; %specify which models to fit diffusivities for -- FH - 1 ; DSM - 2; FH-LM - 3
    
  %specify thermodynamic assumptions
    sysInfo.noThermoCoupling = 0;
    sysInfo.memPhaseModel_FicksOG = 0; %==1 if you want to fit Fick's diffusivities
    
  %specify diffusional relationships **set for parameter fitting**   
    sysInfo.diffModel = "NoCoupling";  
    sysInfo.swlDiffModel = "none";
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%numerical method sepcs **best specs for fitting code**
    sysInfo.numMethod = "NormShootAlg"; %classical shooting algorithm (SA) [Recommended]
  %solver specs
    sysInfo.solverSpec = "trust-region-dogleg"; %default, use for FUD
    sysInfo.iterDetail = 0; %0 = all main solver output will be suppressed 
  %initial guess specs
    sysInfo.initGuessApprox = 0; %0 = use feed mol fraction vector and small total flux value (use if IVP approx fails)
    sysInfo.nodalGuessApprox = 0;
    sysInfo.customInitGuessApprox = 0; %set to match your experimental compostion number plus total flux as [n+1 by 1] vector for custom initial guess
  %sys of eqns spec
    sysInfo.currentStateLit_eqnSetup = 0; %default == 0 (no difference to genreal simulations)
    sysInfo.noGammaFugacityODEs = 0; % ONLY USE FOR FULL MS (FICK == 0 and noThermo == 0)! option keep MS system interms of fugacity gradients and do not make chain rule trick to get gamma matrix
%------------------------------------------------------------------------------------------------------------------------------------%    
