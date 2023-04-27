%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalMSA_ODE_RHS(odeState,algState,solverState,params) %
% Description: Evaluate differncial RHS equations.                          %
% Input:       odeState    - n+1 dimensional vector of membrrane phase      %
%                              volume fractions                             %
%              algState    - n dimensional vector of penetrant fugacities   %
%                              (torr)                                       %
%              solverState - n+1 dimensional variable shooting point values %
%                              and usual localCompFlux values for current   % 
%                              iteration                                    %
%              params      - struct of system parameters                    % 
%                           (see dataBank function for specs)               %                        
% Output:      odeRHS      - function value of ODE RHS for DAE solver       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [odeRHS] = asyMemLocalODE_RHS_casADi(odeState,algState,solverState,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    import casadi.*
    n = params.n;
    odeRHS = SX.sym('G',n+1,1);
    localCompFlux = solverState;
    stateVec = [odeState;algState];
    diffs = params.diffs;
    chis = params.chis;
    HanSolParam = params.HanSolParam;
    Vs = params.Vs;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate RHS ODE
[B,invGam] = matrixEvalGammaB_FH_LM_casADi(stateVec,params,solverState);
    odeRHS(1:n) = -invGam*B*(localCompFlux(1:n)...
        .*Vs(1:n).*localCompFlux(n+1)./sum(Vs(1:n).*localCompFlux(1:n))*1000/3600); %note correct convertsion from LMH to UM3/UM2/s
    odeRHS(n+1) = -sum(odeRHS(1:n));
%------------------------------------------------------------------------------------------------------------------------------------% 
end

    
   
    
    
     
            

    

            
    
            
        
    
        


