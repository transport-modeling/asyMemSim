%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    vol2molFrac_casADi(stateVar,params)                        %
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

function [molFracStateVar] = vol2molFrac_casADi(stateVar,params)
    import casadi.*
    n = params.n;
    molFrac = SX.sym('xFrac',n,1);
    Vs = params.Vs;
    n = params.n;
    phis = stateVar(1:n)/sum(stateVar(1:n));
    funRHS = molFrac-phis.*(sum(Vs(1:n).*molFrac))./Vs(1:n);
    funRHS(n) = 1-sum(molFrac);
    rfp = struct('x', molFrac, 'g', funRHS );
    solver_opts = struct('use_preconditioner',1,'exact_jacobian',1,'abstol',1E-5);
    rf = rootfinder('rf', 'kinsol', rfp, solver_opts);
    molFracStateVar = rf(params.yf, 0);
%------------------------------------------------------------------------------------------------------------------------------------% 
end