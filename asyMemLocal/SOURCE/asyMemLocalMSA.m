%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalMSA(phiFeed,localCompFluxGuess,fugFeed,params)  %
% Description: Initialize multiple shooting point algorithm (MSA) for     % 
%                asyMem local flux.                                       %
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

function [localCompFlux,exitFlag] = asyMemLocalMSA(phiFeed,localCompFluxGuess,fugFeed,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%MSA init
    if params.casADi == 1
      %casADi MSA route
        N = params.numShootPoints; %number of shooting points
        if params.pervapMode == 0
            [localCompFlux] = asyMemLocalMSA_casADi(phiFeed,localCompFluxGuess,fugFeed,params,N);
        else
            [localCompFlux] = asyMemLocalMSA_casADi_pervap(phiFeed,localCompFluxGuess(1:params.n),fugFeed,params,N);
        end
        exitFlag = 'seeAboveSolverOutput';
    else
      %fsolve MSA route   
        N = params.numShootPoints; %number of shooting points
        if params.memPhaseModel == 1
            for i = 1:N-1
                localCompFluxGuess = [localCompFluxGuess;phiFeed*(N-i+1+0.001)/N];
            end
        elseif params.memPhaseModel == 2 || params.memPhaseModel == 3
            for i = 1:N-1
                localCompFluxGuess = [localCompFluxGuess;(phiFeed)*(N-i+1+0.001)/N...
                    ;fugFeed*(N-i+1+0.001)/N];
            end    
        end
        if params.iterDetail == 1
            options = optimoptions(@fsolve,'Display','iter','FunctionTolerance',1E-8,'MaxFunctionEvaluations',45000000,...
                'MaxIterations',6000000,'Algorithm',params.solverSpec);
        else
            options = optimoptions(@fsolve,'Display','off','FunctionTolerance',1E-6,'MaxFunctionEvaluations',45000,...
                'MaxIterations',6000,'Algorithm',params.solverSpec);
        end
        if params.pervapMode == 0
            [localCompFlux,~,exitFlag,~,J] = fsolve(@(localCompFlux)asyMemLocalMSA_RHS(phiFeed,localCompFlux,fugFeed,params,N),...
                localCompFluxGuess,options);
        else
            [localCompFlux,~,exitFlag,~,J] = fsolve(@(localCompFlux)asyMemLocalMSA_RHS_pervap(phiFeed,localCompFlux,fugFeed,params,N),...
                [localCompFluxGuess(1:params.n);localCompFluxGuess(params.n+2:end)],options);
            
    end
%------------------------------------------------------------------------------------------------------------------------------------% 
end

%% fsolve options

%,'OutputFcn',@outfun,
%,'Algorithm','levenberg-marquardt'
%'Algorithm','trust-region-dogleg'
%'Algorithm','trust-region','PrecondBandWidth',params.n+1,'MaxPCGiter',3000,
