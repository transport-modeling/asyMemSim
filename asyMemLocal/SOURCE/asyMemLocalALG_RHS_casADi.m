%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalMSA_ALG_RHS(odeState,algState,params)           %
% Description: Evaluate algebraic RHS equations.                          %
% Input:       odeState - n+1 dimensional vector of membrrane phase       %
%                           volume fractions                              %
%              algState - n dimensional vector of penetrant fugacities    %
%                           (torr)                                        %
%              params   - struct of system parameters                     % 
%                           (see dataBank function for specs)             % 
%              solverState - n+1 dim variable shooting point values       %
%                              and usual localCompFlux values for current % 
%                              iteration                                  %
% Output:      algRHS   - function value of algebraic RHS for DAE solver  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [algRHS] = asyMemLocalALG_RHS_casADi(odeState,algState,params,solverState)
%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    import casadi.*
    n = params.n;
    stateVec = [odeState;algState];
    T = params.T;
    R = params.R;
    psat = params.psat;
    Vs = params.Vs;
    Pu = params.Pu; 
%------------------------------------------------------------------------------------------------------------------------------------% 
    algRHS = evalFH_LM_RHS_casADi(stateVec,params,solverState)-algState./(psat.*exp(Vs(1:n).*(Pu-psat*0.00131)/(R*T)));
%------------------------------------------------------------------------------------------------------------------------------------% 
end
    
   
    
    
     
            

    

            
    
            
        
    
        


