%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    vol2molFrac_RHS(stateVar,molFrac,params)                   %
% Description: Solve for membrane phase volume fractions of feed side of  %
%                membrane based on Flory-Huggins (FH) sorption model      % 
% Input:       stateVar        - (ODE, FH) n+1 dimensional vector of vol  %
%                                  fractions in membrane phase            %
%                                (DAE, FH-LM or DSM) 2*n+1 dimensinal vec %
%                                  of n+1 vol fractions and n fugacities  %
%                                  of membrane phase                      %
%              molFrac         - vector of mol fractions                  %
%              params          - struct of system parameters              %
%                                 (see dataBank function for specs)       %
% Output:      funRHS - function values vector for FH equations           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funRHS] = vol2molFrac_RHS(stateVar,molFrac,params)
    Vs = params.Vs;
    n = params.n;
    funRHS = zeros(n,1);
    phis = stateVar(1:n)/sum(stateVar(1:n));
    for i=1:n-1
    funRHS(i) = molFrac(i)-phis(i)*(sum(Vs(1:n).*molFrac))/Vs(i);
    end
    funRHS(n) = 1-sum(molFrac);
end