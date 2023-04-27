%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalFUD(phiFeed,localCompFluxGuess,fugFeed,params)  %
% Description: Initialize full discritization (FUD) algortihm for solving %
%                asyMem local flux model.                                 %
% Input:       phiFeed            - n+1 dimensional vector of feed side   %
%                                      membrane phase volume fractions    %
%              localCompFluxGuess - guess vector of localCompFlux for     %
%                                     solver                              %
%              fugFeed            - n dimensional vector of penetrant     %
%                                     feed side fugacities (torr)         %
%              params             - struct of system parameters           %
%                                     (see dataBank function for specs)   %
% Output:      localCompFLux      - n+1 dimensional vector of support     %
%                                    layer compositions and total local   %
%                                    mem flux                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [localCompFlux,exitFlag,output] = asyMemLocalFUD(phiFeed,localCompFluxGuess,fugFeed,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%finite differences method for solving MS BVP --- full discritization (FUD)
    n = params.n;
    N = params.numNodes; %number of grid points
    chis = params.chis;
    diffs = params.diffs;
    zspan = [0 params.lmem];
    if params.memPhaseModel == 1 && params.noGammaFugacityODEs == 0
        solVecGuess = zeros((n+1)*(N)+n+1-n,1);
        if params.nodalGuessApprox == 0
            0;
        else
            diffVar = 0;
            funIVP = @(z,stateVar)asyMemLocal_IVP(z,stateVar,diffVar,localCompFluxGuess,params); %solve IVP
            [~,stateVars] = ode15s(funIVP,zspan,phiFeed);
            numIntPoints = size(stateVars,1);
        end
        for i = 1:n+1
            if  N == 1
                solVecGuess(1+i:n+1:(n+1)*(N)+n+1-n) = phiFeed(i);
            else
                if params.nodalGuessApprox == 0
                    if i==n+1
%                         tempVecPhi = linspace(phiFeed(i),0.999,N).';
                        tempVecPhi = linspace(phiFeed(i),phiFeed(i),N).';
                    else
%                         tempVecPhi = linspace(phiFeed(i),0.001*phiFeed(i),N).';
                        tempVecPhi = linspace(phiFeed(i),phiFeed(i),N).';
                    end
                    solVecGuess(1+i:n+1:(n+1)*(N)+n+1-n) = [tempVecPhi];
                else %==1
                    solVecGuess(1+i:n+1:(n+1)*(N)+n+1-n) =  ...
                        interp1(stateVars(:,i),1:(numIntPoints-1)/(N-1):numIntPoints);
                end
            end
        end
        solVecGuess(1) = localCompFluxGuess(n+1);
    elseif params.memPhaseModel == 2 || params.memPhaseModel == 3 || (params.memPhaseModel == 1 && params.noGammaFugacityODEs == 1)
        solVecGuess = zeros((n+1+n)*(N)-n+n+1,1);
        if params.nodalGuessApprox == 0
            0;
        else % == 1
            funIVP = @(z,W,D)asyMemLocal_IVP(z,W,D,localCompFluxGuess,params); %solve IVP
            D0 = ones(n+1+n,1); %FH-DSM
            opt = odeset('InitialSlope', D0,'RelTol',1e-4,'AbsTol',1e-6);
            [W0_new,D0_new] = decic(funIVP,0,[phiFeed;fugFeed],[ones(n+1,1);zeros(n,1)],D0,zeros(1,n+1+n),opt); %FH-DSM
            opt = odeset(opt,'InitialSlope', D0_new,'RelTol',1e-4,'AbsTol',1e-6);
            [~,stateVars] = ode15i(funIVP,zspan,W0_new,D0_new,opt);
            numIntPoints = size(stateVars,1);
        end
        for i = 1:n+1
            if i == n+1 && params.nodalGuessApprox == 0
                tempVecPhi = linspace(phiFeed(i),0.999,N).';
            elseif i ~= n+1
                if N == 1
                    solVecGuess(n+1+n+1+i:n+1+n:(n+1+n)*(N)-n+n+1) = fugFeed(i);
                else
                    if params.nodalGuessApprox == 0
                        tempVecPhi = linspace(phiFeed(i),0.001*phiFeed(i),N).';
                        tempVecFug = linspace(fugFeed(i),0.001*fugFeed(i),N-1).';
                        solVecGuess(n+1+n+1+i:n+1+n:(n+1+n)*(N)-n+n+1) = [tempVecFug];
                    else %==1
                        solVecGuess(n+1+n+1+i:n+1+n:(n+1+n)*(N)-n+n+1) =...
                            interp1(stateVars(:,i+n+1),1:(numIntPoints-2)/(N-2):numIntPoints-1);
                    end
                end
            end
            if N == 1
                solVecGuess(n+1+i:n+1+n:(n+n+1)*(N)+n+1) = phiFeed(i);
            else
                if params.nodalGuessApprox == 0
                    solVecGuess(n+1+i:n+1+n:(n+n+1)*(N)+n+1) = [tempVecPhi];
                else %==1
                    solVecGuess(n+1+i:n+1+n:(n+n+1)*(N)+n+1) = ...
                        interp1(stateVars(:,i),1:(numIntPoints-1)/(N-1):numIntPoints);
                end
            end
        end 
        solVecGuess(1:n+1) = localCompFluxGuess;
    end
    funRHS = @(sol_vec)asyMemLocalFUD_RHS(phiFeed,fugFeed,sol_vec,N,params);
    if params.iterDetail == 1
        options = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1E-8,'MaxFunctionEvaluations',45000000,...
            'MaxIterations',6000000,'Algorithm',params.solverSpec);
    else
        options = optimoptions(@fsolve,'Display','off','FunctionTolerance',1E-8,'MaxFunctionEvaluations',45000,...
            'MaxIterations',6000,'Algorithm',params.solverSpec);
    end
    [solVec,~,exitFlag,output,J] = fsolve(funRHS,solVecGuess,options);
    if params.memPhaseModel == 1 && params.noGammaFugacityODEs == 0
        phiFinal = solVec((n+1)*N-n+1:(n+1)*N+n+1-n);
        [params.chis,~] = correlationEval(phiFinal,diffs,chis,params);
        localCompFlux(1:n) = phis2yPhaseEq_FH_RHS(phiFinal,params);
        localCompFlux(n+1) = solVec(1);
        localCompFlux = localCompFlux.';
    elseif params.memPhaseModel == 2 || params.memPhaseModel == 3 || (params.memPhaseModel == 1 && params.noGammaFugacityODEs == 1)
        localCompFlux = solVec(1:n+1);
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

end
%,'OutputFcn',@outfun,