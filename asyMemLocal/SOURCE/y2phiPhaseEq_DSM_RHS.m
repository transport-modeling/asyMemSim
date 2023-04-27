%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    y2phiPhaseEq_FH_RHS(phis,params)                           %
% Description: Evaluate RHS function for solving of feed side membrane    %
%              phase vol frac based on Flory-Huggins (FH) sorption model. % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
%              phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
% Output:      funRHS - function values vector for FH equations           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funRHS] = y2phiPhaseEq_DSM_RHS(phis,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    ys = params.yf;
    n = params.n;
    Vs = params.Vs;
    psat = params.psat;
    R = params.R;
    T = params.T;
    Pu = params.Pu;
    ks = params.ks;
    bs = params.bs;
    Ch = params.Ch;
    fugFeed_IS = ys.*psat.*exp(Vs(1:n).*(Pu-psat*0.00131)/(R*T));
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate RHS of DSM
    for i = 1:n
        funRHS(i) = -phis(i)+(ks(i)*fugFeed_IS(i) + Ch(i)*bs(i)*fugFeed_IS(i)/(1+sum(bs.*fugFeed_IS)))*phis(n+1);
    end
    funRHS(n+1) = 1-sum(phis(1:n+1));
%------------------------------------------------------------------------------------------------------------------------------------% 

end
   
    