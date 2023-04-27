%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    phi2fugPhaseEq_FH_LM(params,phis)                          %
% Description: Solve for membrane phase fugacities based on combined      %
%                Flory-Huggins and Langmuir volume fractions              %
%                (FH-LM) sorption model                                   % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
%              phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
% Output:      fug    - n dimensional vector of membrane phase fugacities %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [fug] = phi2fugPhaseEq_FH_LM(params,phis)

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve for phi based on bulk feed composition
    n = params.n;
    fugGuess(1:n) = params.yf.*0.001;
    options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',50000);
    fug = fsolve(@(fug)DAEevalFH_LM_RHS([phis;fug],params),fugGuess.',options);
%------------------------------------------------------------------------------------------------------------------------------------% 

end

