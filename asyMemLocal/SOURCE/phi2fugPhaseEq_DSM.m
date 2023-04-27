%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    phi2fugPhaseEq_DSM(params,phis)                            %
% Description: Solve for membrane phase fugacities based on Dual-Sorption %
%                Mode (DSM) sorp model volume fractions                   % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
%              phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
% Output:      fug    - n dimensional vector of membrane phase fugacities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fug] = phi2fugPhaseEq_DSM(params,phis)

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve for phi based on bulk feed composition
    n = params.n;
    fugGuess(1:n) = params.yf.*0.001;
    options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000);
    fug = fsolve(@(fug)DAEevalDSM_RHS([phis;fug],params),fugGuess.',options);
%------------------------------------------------------------------------------------------------------------------------------------% 

end

