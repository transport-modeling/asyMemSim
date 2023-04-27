function [sysInfo] = solverModelSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
    sysInfo.numMethod = "NormShootAlg"; %classical shooting algorithm (SA) [Recommended]
  %solver specs
    sysInfo.solverSpec = "trust-region-dogleg"; %default, use for FUD
    sysInfo.iterDetail = 0; %0 = all main solver output will be suppressed 
  %initial guess specs
    sysInfo.initGuessApprox = 1; %0 = use feed mol fraction vector and small total flux value (use if Forward Euler IVP approx fails)
    sysInfo.nodalGuessApprox = 0;
    sysInfo.customInitGuessApprox = []; %set to match your experimental compostion number plus total flux as [n+1 by 1] vector for custom initial guess
  %sys of eqns spec
    sysInfo.currentStateLit_eqnSetup = 0; %default == 0 (no difference to genreal simulations)
    sysInfo.noGammaFugacityODEs = 0; % ONLY USE FOR FULL MS and SA or FUD method (FICK == 0 and noThermo == 0)! option keep MS system interms of fugacity gradients and do not make chain rule trick to get gamma matrix
%------------------------------------------------------------------------------------------------------------------------------------%    
