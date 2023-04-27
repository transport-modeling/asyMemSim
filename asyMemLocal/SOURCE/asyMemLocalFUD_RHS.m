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

function [funRHS] = asyMemLocalFUD_RHS(phiFeed,fugFeed,solVec,N,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    n = params.n;
    diffVar = zeros(n+1+n,1);
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
    if N ~= 1
        z = linspace(0,lmem,(N)).'; %genreate mesh
        h = z(2)-z(1);
    else
        z = lmem;
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate RHS for discritized BVP
    if params.memPhaseModel == 1 && params.noGammaFugacityODEs == 0
        funRHS = zeros((n+1)*N+n+1-n,1);
        phiFinal = solVec((n+1)*N-n+1:(n+1)*N+n+1-n);
        [params.chis,~] = correlationEval(phiFinal,diffs,chis,params);
        localCompFlux(1:n) = phis2yPhaseEq_FH_RHS(phiFinal,params);
        localCompFlux(n+1) = solVec(1);
        localCompFlux = localCompFlux.';
        funRHS(1) = 1-sum(localCompFlux(1:n));
        for i = 2:n+1:(n+1)*N+n+1-n  %build nodal equation space for FUD method
            nodeStatePhi_i = solVec(i:i+n);
            [ODE_RHS_i] = asyMemLocal_IVP(z,nodeStatePhi_i,diffVar,localCompFlux,params);
            if i == 2
                %Interior BC 
                   funRHS(i:i+n) = nodeStatePhi_i-phiFeed;
                %Krishnas (and Izak) midpoint method (N = 1)
                if N == 1
                    funRHS(i:i+n) = -(solVec(i:i+n)-phiFeed)+...
                        lmem*(asyMemLocal_IVP(z,(nodeStatePhi_i+phiFeed)/2,diffVar,localCompFlux,params));
                end
            elseif i == (n+1)*N-n+n+1-n
                %BD
                funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1:i-1))+...
                    h*ODE_RHS_i;
            else
                %CD
               funRHS(i:i+n) = -(solVec(i+n+1:i+2*n+1)-solVec(i-n-1:i-1))+ ...
                   2*h*ODE_RHS_i;
            end    
        end
    elseif params.memPhaseModel == 2 || params.memPhaseModel == 3 || (params.memPhaseModel == 1 && params.noGammaFugacityODEs == 1)
        localCompFlux = solVec(1:n+1);
        funRHS = zeros((n+1+n)*N+n+1-n,1);
        for i = n+2:n+1+n:(n+n+1)*N+n+1-n
            nodeStatePhi_i = solVec(i:i+n); 
            if N == 1 || i == (n+n+1)*N+n+1-n-n
                nodeStateFug_i = localCompFlux(1:n).*(exp(-Vs(1:n).*(Pu-Pd)/(R*T))).*params.purFug;
            else
                nodeStateFug_i = solVec(i+(n+1):i+n+n);
            end
            [RHS] = asyMemLocal_IVP(z,[nodeStatePhi_i;nodeStateFug_i]...
                ,diffVar,localCompFlux,params);
            ODE_RHS_i = RHS(1:n+1);
            ALG_RHS_i = RHS(n+2:end);
            if i == n+2
                %Interior BC 
                if N ~= 1
                   funRHS(i:i+n+n) = [nodeStatePhi_i;nodeStateFug_i]-[phiFeed;fugFeed];
%                    funRHS(i:i+n+n) = [nodeStatePhi_i;exp(nodeStateFug_i)]-[phiFeed;fugFeed];
                end

                %Krishnas (and Izak) midpoint method (N = 1)
                if N == 1
                    [RHS_N_1] = asyMemLocal_IVP(z,[(nodeStatePhi_i+phiFeed)/2;(nodeStateFug_i+fugFeed)/2]...
                        ,diffVar,localCompFlux,params);
                    ODE_RHS_i_N_1 = RHS_N_1(1:n+1);
                    ALG_RHS_i_N_1 = RHS_N_1(n+2:end);
                    if params.noGammaFugacityODEs == 0
                        funRHS(i:i+n) = -(nodeStatePhi_i-phiFeed)+lmem*ODE_RHS_i_N_1;
                    else
                        funRHS(i:i+n) = -(nodeStatePhi_i+phiFeed)/2.*[(log(nodeStateFug_i)-log(fugFeed));0]+lmem*ODE_RHS_i_N_1;
                        funRHS(i+n) = funRHS(i+n)+1-sum(nodeStatePhi_i);
                    end
                end
            elseif i == (n+n+1)*N+n+1-n-n
                if params.noGammaFugacityODEs == 1
                %BD
                    funRHS(i:i+n) = -solVec(i:i+n).*[log(nodeStateFug_i)-log(solVec(i-n-1-n+n+1:i-1-n+n));0]+...
                        h*ODE_RHS_i;
                    funRHS(i+n) = funRHS(i+n)+1-sum(solVec(i:i+n));
                %BD
                else
                    funRHS(i:i+n) = -(solVec(i:i+n)-solVec(i-n-1-n:i-1-n))+...
                        h*ODE_RHS_i;
                    funRHS(i+n) = funRHS(i+n)-sum(ODE_RHS_i(1:n))*h;
                end
            else
                if params.noGammaFugacityODEs == 1
                %CD
                    if i == (n+n+1)*N+n+1-n-n-(n+1+n)
                    funRHS(i:i+n) = -solVec(i:i+n).*[(log(localCompFlux(1:n).*(exp(-Vs(1:n).*(Pu-Pd)/(R*T))).*params.purFug)-log(solVec(i-n-1-n+n+1:i-1-n+n)));0]+ ...
                        2*h*ODE_RHS_i;
                    else
                        funRHS(i:i+n) = -solVec(i:i+n).*[(log(solVec(i+n+1+n+n+1:i+2*n+1+n+n))-log(solVec(i-n-1-n+n+1:i-1-n+n)));0]+ ...
                            2*h*ODE_RHS_i;
                    end
                    funRHS(i+n) = funRHS(i+n)+1-sum(solVec(i:i+n));
                else
                %CD
                    funRHS(i:i+n) = -(solVec(i+n+1+n:i+2*n+1+n)-solVec(i-n-1-n:i-1-n))+ ...
                        2*h*ODE_RHS_i;
                    funRHS(i+n) = funRHS(i+n)-sum(ODE_RHS_i(1:n))*2*h;
                end
            end 
            if N == 1 || i == (n+n+1)*N+n+1-n-n 
                [RHS_N_1] = asyMemLocal_IVP(z,[nodeStatePhi_i;nodeStateFug_i]...
                        ,diffVar,localCompFlux,params);
                funRHS(1:n) = RHS_N_1(n+2:end);
            elseif i ~= n+2
                funRHS(i+n+1:i+n+n) = ALG_RHS_i;
            end
        end
    %summation equation for molar composition in support layer
    funRHS(n+1) = 1-sum(localCompFlux(1:n));
    end

%------------------------------------------------------------------------------------------------------------------------------------% 
end




