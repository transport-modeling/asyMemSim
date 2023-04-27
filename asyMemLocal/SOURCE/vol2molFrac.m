%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    vol2molFrac(stateVar,params)                               %
% Description: Solve for membrane phase volume fractions of feed side of  %
%                membrane based on Flory-Huggins (FH) sorption model      % 
% Input:       stateVar        - (ODE, FH) n+1 dimensional vector of vol  %
%                                  fractions in membrane phase            %
%                                (DAE, FH-LM or DSM) 2*n+1 dimensinal vec %
%                                  of n+1 vol fractions and n fugacities  %
%                                  of membrane phase                      %
%              params          - struct of system parameters              %
%                                 (see dataBank function for specs)       %
% Output:      molFracStateVar - vector of mol fractions                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [molFracStateVar] = vol2molFrac(stateVar,params)
    options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000);
    molFracStateVar = fsolve(@(molFrac)vol2molFrac_RHS(stateVar,molFrac,params),zeros(params.n,1),options);
end