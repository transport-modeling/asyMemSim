%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalFUD_RHS(phiFeed,solVec,N,params)                %
% Description: Evaluate RHS function for FUD solver.                      %
% Input:       phiFeed - n+1 dimensional vector of feed side              %
%                          membrane phase volume fractions                %
%              fugFeed - n dimensional vector of feed side fugacities     %
%              solVec  - vector of nodal values and localCompFlux for     %
%                          current iteration                              %
%              N       - number of discritized nodes                      %
%              params  - struct of system parameters                      %
%                          (see dataBank function for specs)              %
% Output:      funRHS  - function value vector for nonlinear solver       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funRHS] = asyMemLocalFUD_RHS_pervap(phiFeed,fugFeed,solVec,N,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    n = params.n;
    diffVar = zeros(n+1,1);
    if params.memPhaseModel == 1 || params.memPhaseModel == 3
        chis = params.chis; %note chi_ji = chi_ij
    elseif params.memPhaseModel == 2
        chis = zeros(params.n+1);
    end
    Pu = params.Pu;
    Pd = params.Pd;
    Vs = params.Vs;
    R = params.R;
    T = params.T;
    diffs = params.diffs;
    psat = params.psat;
    HanSolParam = params.HanSolParam;
    lmem = params.lmem;
    z = linspace(0,lmem,(N)).'; %genreate mesh
    h = z(2)-z(1); 
    localCompFlux = solVec(1:n+1);
    %FH
    funRHS = zeros((n+1)*N+n+1,1);
    %FH-LM or DSM
    %funRHS = zeros((n+1+n)*N+n+1,1);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate RHS for discritized BVP
    if params.memPhaseModel == 1
        for i = n+2:n+1:(n+1)*N+n+1  %build nodal equation space for FUD method
            nodeStatePhi_i = solVec(i:i+n);
            [ODE_RHS_i] = asyMemLocal_IVP(z,nodeStatePhi_i,diffVar,localCompFlux,params);
            if i == n+2
                %Interior BC 
                   funRHS(i:i+n) = nodeStatePhi_i-phiFeed;
                %CD
%                  funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-phiFeed)+ ...
%                      2*h*ODE_RHS_i;
                %BD
%                 funRHS(i:i+n) = -(solVec(i:i+n)-phiFeed)+...
%                    h*ODE_RHS_i;
                %Krishnas (and Izak) midpoint method (N = 1)
%                 funRHS(i:i+n) = -(solVec(i:i+n)-phiFeed)+...
%                    lmem*(asyMemLocal_IVP(z,(nodeStatePhi_i+phiFeed)/2,diffVar,localCompFlux,params));
            elseif i == (n+1)*N-n+n+1
                %BD
                funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1:i-1))+...
                    h*ODE_RHS_i;
            else
                %CD
               funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-solVec(i-n-1:i-1))+ ...
                   2*h*ODE_RHS_i;
                %BD
    %              funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1:i-1))+...
    %                  h*ODE_RHS_i;
                %FD
%                 funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-solVec(i:i+n))+...
%                    h*ODE_RHS_i;
            end    
        end
        phiFinal = solVec((n+1)*N-n+n+1:(n+1)*N+n+1);
        [params.chis,~] = correlationEval(phiFinal,diffs,chis,params);
        ypFinal = phis2yPhaseEq_FH_RHS(phiFinal,params);
    elseif params.memPhaseModel == 2 || params.memPhaseModel == 3
        for i = n+2:n+1+n:(n+n+1)*N+n+1
            nodeStatePhi_i = solVec(i:i+n); 
            nodeStateFug_i = solVec(i+(n+1):i+n+n);
            [RHS] = asyMemLocal_IVP(z,[nodeStatePhi_i;nodeStateFug_i]...
                ,diffVar,localCompFlux,params);
            ODE_RHS_i = RHS(1:n+1);
            ALG_RHS_i = RHS(n+2:end);
            if i == n+2
                %Interior BC 
                   funRHS(i:i+n+n) = [nodeStatePhi_i;nodeStateFug_i]-[phiFeed;fugFeed];
                %CD
    %              funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-phiFeed)+ ...
    %                  2*h*ODE_RHS_i;
                %BD
    %             funRHS(i:i+n) = -(solVec(i:i+n)-phifs)+...
    %                h*ODE_RHS_i;
                %Krishnas (and Izak) midpoint method (N = 1)
    %             funRHS(i:i+n) = -(solVec(i:i+n)-phifs)+...
    %                lmem*(asyMemLocal_IVP(z,[(nodeStatePhi_i+phiFeed)/2;...
    %                (nodeStateFug_i+fugFeed)/2],diffVar,localCompFlux,params));
            elseif i == (n+n+1)*N+n+1-n-n
                %BD
                funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1-n:i-1-n))+...
                    h*ODE_RHS_i;
                funRHS(i+n) = funRHS(i+n)-sum(ODE_RHS_i(1:n))*h;
            else
                %CD
               funRHS(i:i+n) = -(solVec(i+n+1+n:i+2*n+1+n)-solVec(i-n-1-n:i-1-n))+ ...
                   2*h*ODE_RHS_i;
                %BD
    %              funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1:i-1))+...
    %                  h*ODE_RHS_i;
                %FD
%                 funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-solVec(i:i+n))+...
%                    h*ODE_RHS_i;
                funRHS(i+n) = funRHS(i+n)-sum(ODE_RHS_i(1:n))*2*h;
            end 

            funRHS(i+n+1:i+n+n) = ALG_RHS_i;
        end
        fugFinal = solVec((2*n+1)*N-n+1+n+1:(2*n+1)*N-n+n+1+n);
        ypFinal = fugFinal./(psat.*exp(-Vs(1:n).*(Pu-2*Pd+psat*0.00131)/(R*T)));
    end
    %n+1 eqns from previous methods
    funRHS(1:n) = localCompFlux(1:n)-ypFinal;
    funRHS(n+1) = 1-sum(ypFinal);
%------------------------------------------------------------------------------------------------------------------------------------% 

end

