%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalFoEu_QUICK(phiFeed,localCompFluxGuess,fugFeed   %
%                ,params)                                                 %
% Description: Initialize shooting algorithm for QUICK asyMem local flux  %
%                solve for quick quess using forward IVP Euler's method.  %
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

function [localCompFlux] = asyMemLocalFoEu_QUICK(phiFeed,localCompFluxGuess,fugFeed,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%MATLAB fsolve route     
    options = optimoptions(@fsolve,'Display','off','FunctionTolerance',1E-8,'MaxFunctionEvaluations',...
        4000000,'MaxIterations',50000,'OptimalityTolerance',1E-8);
    [localCompFlux,~,~,~,J] = fsolve(@(localCompFlux)asyMemLocalFoEu_QUICK_RHS(phiFeed,localCompFlux,fugFeed,params),...
        localCompFluxGuess,options);
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