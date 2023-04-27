%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    evalDSM_RHS(stateVar,params)                               %
% Description: Evaluate RHS of Dual-mode sorption equations.              %
% Input:       stateVar - (ODE, FH) n+1 dimensional vector of volume      %
%                           fractions in membrane phase                   %
%                         (DAE, FH-LM or DSM) 2*n+1 dimensinal vector     %
%                           of n+1 volume fractions and n fugacities      %
%                           of membrane phase                             %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
% Output:      DSM_RHS  - nonlinear function value of DSM RHS equation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [DSM_RHS] = DAEevalDSM_RHS(stateVar,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
n = params.n;
bs = params.bs;
ks = params.ks;
Ch = params.Ch;
fs = stateVar(n+2:end);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate RHS
DSM_RHS = stateVar(n+1)*((ks.*fs) + Ch.*bs.*fs/(1+sum(bs.*fs)))-stateVar(1:n);
%------------------------------------------------------------------------------------------------------------------------------------% 

end

