%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    MSAevalFH_RHS(stateVar,params)                             %
% Description: Evaluate RHS of Flory-Huggings equation for MSA.           %
% Input:       stateVar - (ODE, FH) n+1 dimensional vector of volume      %
%                           fractions in membrane phase                   %
%                         (DAE, FH-LM or DSM) 2*n+1 dimensinal vector     %
%                           of n+1 volume fractions and n fugacities      %
%                           of membrane phase                             %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
%              solverState - n+1 dim variable shooting point values       %
%                              and usual localCompFlux values for current % 
%                              iteration                                  %
% Output:      FH_RHS   - nonlinear function value of FH RHS equation     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FH_LM_RHS] = evalFH_LM_RHS_casADi(stateVar,params,solverState)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parmeters and intialize/solve for variables
    import casadi.*
    n = params.n;   
    phiFH = SX.sym('pFH',n+1,1);
    phiFHfun = SX.sym('fFH',n,1);
    Ch = params.Ch;
    bs = params.bs;
    Vs = params.Vs;
    fs = stateVar(n+2:end);
    phiFH(1:n,1) = (stateVar(1:n)-stateVar(n+1).*Ch.*bs.*fs/(1+sum(bs.*fs)));
    phiFH(n+1,1) = 1-sum(phiFH(1:n));
    sxDiffs = SX.sym('sxDiffs',n,n+1);
    sxChis = SX.sym('sxChis',n+1,n+1);
    sxChis(n+1,:) = params.chis(n+1,:);
    sxChis(:,n+1) = params.chis(:,n+1);
    if n==1
        sxDiffs(1,n+1) = params.diffs(1,n+1);
    else
        sxDiffs(:,n+1) = params.diffs(1:n,n+1);
    end
    [sxChis,~] = correlationEval_casADi(stateVar,sxDiffs,sxChis,params,solverState);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evalute FH RHS
    for i = 1:n
        phiFHfun(i) = (1-phiFH(i))-Vs(i)*(sum(phiFH./Vs)-phiFH(i)/Vs(i))+Vs(i)*(sum(phiFH(1:i-1)...
            .*(sxChis(1:i-1,i))./Vs(1:i-1)))*(sum(phiFH)-phiFH(i))+sum(phiFH(1+i:n+1).*((sxChis(i,1+i:n+1)).'))*(sum(phiFH)-phiFH(i));        
        for j = 1:n
            if j == i, continue, end
            phiFHfun(i) = phiFHfun(i)-Vs(i)/Vs(j)*sum(((sxChis(j,j+1:n+1)).').*phiFH(j+1:n+1))*phiFH(j);
        end  
    end
    FH_LM_RHS = phiFH(1:n).*exp(phiFHfun);
%------------------------------------------------------------------------------------------------------------------------------------% 

end

