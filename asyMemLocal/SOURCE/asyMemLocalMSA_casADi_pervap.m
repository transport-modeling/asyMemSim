%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalMSA(phiFeed,localCompFluxGuess,fugFeed,params)  %
% Description: Initialize multiple shooting point algorithm (MSA) for     % 
%                asyMem local flux.                                       %
% Input:       phiFeed            - n+1 dimensional vector of feed side   %
%                                      membrane phase volume fractions    %
%              localCompFluxGuess - guess vector of localCompFlux for     %
%                                     solver                              %
%              fugFeed            - n dimensional vector of penetrant     %
%                                     feed side fugacities (torr)         %
%              params             - struct of system parameters           %
%                                     (see dataBank function for specs)   %
% Output:      localCompFLux      - n+1 dimensional vector of support     %
%                                    layer compositions and total local   %
%                                    mem flux                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [localCompFlux] = asyMemLocalMSA_casADi_pervap(phiFeed,localCompFluxGuess,fugFeed,params,N)
%------------------------------------------------------------------------------------------------------------------------------------% 
%initialize algorithm specs
   %casADi::KINSOL or NLPSOL
     import casadi.*
     n = params.n; %number of permeants
%------------------------------------------------------------------------------------------------------------------------------------%      
     
%------------------------------------------------------------------------------------------------------------------------------------%   
%init DAE solve
     odeState = SX.sym('oS',params.n+1,1);
     algState = SX.sym('aS',params.n,1);
     solverState = SX.sym('sS',params.n,1);
     odeRHS = asyMemLocalODE_RHS_casADi_pervap(odeState,algState,solverState,params);
     algRHS = asyMemLocalALG_RHS_casADi(odeState,algState,params,solverState);   
     dae = struct('x', odeState, 'z', algState,'p',solverState, 'ode', odeRHS, 'alg', algRHS);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%init outer solver and integrator instances
     params.solVar = MX.sym('sV',params.n+2*n*(N-1),1);
     for i = 1:N
         if i==1
             opts_init = struct('t0',0,'tf', params.lmem/N,'use_preconditioner',1,'abstol',1E-10,'max_num_steps',20000); 
             inteObj_init = integrator('inteObj', 'idas', dae, opts_init);
             daeSolve_init = inteObj_init('x0',phiFeed,'z0',fugFeed,'p',params.solVar(1:n));
             daeSolveArray = [daeSolve_init];
         else
             opts_i = struct('t0',(i-1)*params.lmem/(N),'tf', (i)*params.lmem/(N),...
                 'use_preconditioner',1,'abstol',1E-10,'max_num_steps',20000); 
             inteObj_i = integrator('inteObj', 'idas', dae, opts_i);
             daeSolve_i = inteObj_i('x0',[params.solVar(n+1+(i-2)*2*n:n+n+(i-2)*2*n);...
                 1-sum(params.solVar(n+1+(i-2)*2*n:n+n+(i-2)*2*n))],'z0',...
                 params.solVar(n+1+n+(i-2)*2*n:n+2*n+(i-2)*2*n),'p',params.solVar(1:n));
             daeSolveArray = [daeSolveArray;daeSolve_i];
         end
     end
     solvRHS = 0.00001 - daeSolveArray(end).zf;
     solvRHS_MultShoot = [];
     for j = 1:N-1
         solvRHS_MultShoot_i = [daeSolveArray(j).xf(1:n)-...
             params.solVar(n+1+(j-1)*2*n:n+n+(j-1)*2*n);...
             daeSolveArray(j).zf-params.solVar(n+1+n+(j-1)*2*n:n+2*n+(j-1)*2*n)];
         solvRHS_MultShoot = [solvRHS_MultShoot;solvRHS_MultShoot_i];
     end
     rfp = struct('x', params.solVar, 'g', [solvRHS;solvRHS_MultShoot]); 
     solver_opts = struct('nlpsol','ipopt');
         %'nlpsol','ipopt'
         %'use_preconditioner',1,'exact_jacobian',1,'abstol',1E-5
         %'constraints',ones(n+1+n,1)
     %solve shooting algorithm
     rf = rootfinder('rf', 'nlpsol', rfp, solver_opts);
     for k = 1:N-1
     localCompFluxGuess = [localCompFluxGuess;phiFeed(1:n);fugFeed];
     end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve MSA
     localCompFlux = rf(localCompFluxGuess, 0);
%------------------------------------------------------------------------------------------------------------------------------------% 
end