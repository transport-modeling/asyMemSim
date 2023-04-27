%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA(phiFeed,localCompFluxGuess,fugFeed,params)   %
% Description: Initialize shooting algorithm for asyMem local flux.       %
% Input:       phiFeed            - n+1 dimensional vector of feed side   %
%                                     membrane phase volume fractions     %
%              localCompFluxGuess - guess vector of localCompFlux for     %
%                                     solver                              %
%              fugFeed            - n dimensional vector of penetrant     %
%                                     feed side fugacities (torr)         %
%              params             - struct of system parameters           %
%                                     (see dataBank function for specs)   %
% Output:      localCompFlux      - n+1 dimensional vector of support     %
%                                     layer compositions and total local  %
%                                     mem flux                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [localCompFlux,exitFlag,output] = asyMemLocalSA(phiFeed,localCompFluxGuess,fugFeed,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%MATLAB fsolve route     
    if params.iterDetail == 1
        options = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1E-8,'MaxFunctionEvaluations',450000,...
            'MaxIterations',600,'Algorithm',params.solverSpec,'OptimalityTolerance',1E-8);
    else
        options = optimoptions(@fsolve,'Display','off','FunctionTolerance',1E-8,'MaxFunctionEvaluations',45000,...
            'MaxIterations',500,'Algorithm',params.solverSpec);
    end
    if params.pervapMode == 0 && params.currentStateLit_eqnSetup == 1
        localCompFluxGuess = localCompFluxGuess(1:params.n)*localCompFluxGuess(params.n+1);
        [localCompFlux,~,exitFlag,output,J] = fsolve(@(localCompFlux)asyMemLocalSA_RHS_OGeqn(phiFeed,localCompFlux,fugFeed,params),...
            localCompFluxGuess,options);
    elseif params.pervapMode == 0 && params.currentStateLit_eqnSetup == 0
        [localCompFlux,~,exitFlag,output,J] = fsolve(@(localCompFlux)asyMemLocalSA_RHS(phiFeed,localCompFlux,fugFeed,params),...
            localCompFluxGuess,options);
    elseif params.pervapMode == 1
        localCompFluxGuess = localCompFluxGuess(1:params.n)*localCompFluxGuess(params.n+1);
        [localCompFlux,~,exitFlag,output,J] = fsolve(@(localCompFlux)asyMemLocalSA_RHS_pervap(phiFeed,localCompFlux,fugFeed,params),...
            localCompFluxGuess,options);
    end
%------------------------------------------------------------------------------------------------------------------------------------%  

end

%% fsolve options

%,'OutputFcn',@outfun,
%,'Algorithm','levenberg-marquardt'
%'Algorithm','trust-region-dogleg'
%'Algorithm','trust-region','PrecondBandWidth',params.n+1,'MaxPCGiter',3000,

%% ect.

   %casADi::KINSOL or NLPSOL
%      import casadi.*
%      n = params.n;
%      
%      %init DAE solve
%      odeState = SX.sym('oS',params.n+1,1);
%      algState = SX.sym('aS',params.n,1);
%      solverState = SX.sym('sS',params.n+1,1);
%      odeRHS = MS_asym_mem_ODE_RHS_casADi(odeState,algState,solverState,params);
%      algRHS = MS_asym_mem_ALG_RHS_casADi(odeState,algState,params);   
%      dae = struct('x', odeState, 'z', algState,'p',solverState, 'ode', odeRHS, 'alg', algRHS);
%      opts = struct('tf', params.lmem,'use_preconditioner',1,'abstol',1E-8,'max_num_steps',20000,'newton_scheme','tfqmr');
%      inteObj = integrator('inteObj', 'idas', dae, opts);
%      
%      %init outer solver
%      params.solVar = MX.sym('sV',params.n+1,1);
%      daeSolve = inteObj('x0',phi_fs,'z0',fug_fs,'p',params.solVar);
%      y_p_final = daeSolve.zf./...
%          (params.p_sat.*exp(-params.V_s(1:n).*(params.P_u-2*params.P_d+params.p_sat)/(params.R*params.T)));
%      solvRHS = params.solVar(1:params.n)-y_p_final;
%      solvRHS_end = 1-sum(y_p_final);
%      rfp = struct('x', params.solVar, 'g', [solvRHS;solvRHS_end]); 
%      solver_opts = struct('nlpsol','ipopt');
%          %'nlpsol','ipopt'
%          %'use_preconditioner',1,'exact_jacobian',1,'abstol',1E-5
%      
%      %solve shooting algorithm
%      rf = rootfinder('rf', 'nlpsol', rfp, solver_opts);
%      comp_perm_flux = rf(comp_perm_flux_guess, 0);
 
   
  %edit for multiple shooting w/ timestepping IVP
%      params.N = 2; %num shooting points
%      for k = 1:params.N-1
%          comp_perm_flux_guess = [comp_perm_flux_guess;phi_fs(1:params.n)];
%      end