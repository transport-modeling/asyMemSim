%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    y2phiPhaseEq_FH_LM(params)                                 %
% Description: Solve for membrane phase volume fractions of feed side of  %
%                membrane based on combined Flory-Huggins and Langmuir    %
%                (FH-LM) sorption model                                   % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
% Output:      phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [phis] = y2phiPhaseEq_FH_LM(params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve for phi based on bulk feed composition
    n = params.n;
    phiGuess(1:n) = params.yf.*0.001;
    phiGuess(n+1) = 1-sum(phiGuess(1:n));
    fugFeed_IS = params.yf.*params.purFug;
    options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000);
    phis = fsolve(@(phi)y2phiPhaseEq_FH_LM_RHS([phi;fugFeed_IS],params),phiGuess.',options);
%------------------------------------------------------------------------------------------------------------------------------------% 



end

