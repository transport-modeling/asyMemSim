%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA_IVP_TimeStep(phiFeed,localCompFlux,params)   %
% Description: Solve initial value problem (IVP) using explicit or        %
%                implicit Euler's method.                                 %
% Input:       phiFeed       - n+1 dimensional vector of feed side        %
%                                      membrane phase volume fractions    %
%              localCompFlux - n+1 dimensional vector of support layer    %
%                                compositions and total local mem flux    %
%              params        - struct of system parameters                %
%                                     (see dataBank function for specs)   %
% Output:      z             - N vector of spacital interval values (um)  %
%              phis          - N x n matrix of membrane phase volume      %
%                                fractions for all values of z            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [z,phis,fugFinal] = asyMemLocalSA_IVP_TimeStep(phiFeed,fugFeed,localCompFlux,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%Euler's Method of Solving Discritized IVP
    n = params.n;
    diffVar = zeros(n+1+n,1);
    N = 5; %number of time points
    z = linspace(0,params.lmem,(N)).'; %genreate mesh
    h = z(2)-z(1); %step size
    phis = zeros(N,params.n+1);
    fugs = zeros(N,params.n);
    phis(1,:) = phiFeed.';
    fugs(1,:) = fugFeed.';
    params.noGammaFugacityODEs = 0;
    for i = 2:N
        prevPhi = phis(i-1,:).';
        if params.memPhaseModel == 1
            prevState = [prevPhi];
        elseif params.memPhaseModel == 2 || params.memPhaseModel == 3
            prevFug = (fugs(i-1,:).')*(N-i-2)/N; %need implicit if using more than 2 steps and FH-LM/DSM
            prevState = [prevPhi;prevFug];            
        end
        %Forward Euler
        dphidz_alg = asyMemLocal_IVP(z,prevState,diffVar,localCompFlux,params);
        phis(i,:) = (prevPhi+dphidz_alg(1:n+1).*h).';
        %Implicit Time-stepping
%         options = optimoptions(@fsolve,'Display','off','Algorithm','levenberg-marquardt');
%         fun = @(currentPhiasyMemLocalSA_IVP_TimeStep_RHS(z,currentPhi,prevPhi,diffVar,localCompFlux,params,h);
% %          phis(i,:) = fsolve(fun,[params.y_f.*.0001;1-sum(params.y_f.*.0001)],options).';
%        phis(i,:) = fsolve(fun,phis(i,:).',options).'; %predictor-corrector
    end
    phiFinal = phis(end,:).';
    if params.memPhaseModel == 2
        fugFinal = phi2fugPhaseEq_DSM(params,phiFinal);
    elseif params.memPhaseModel == 3
        fugFinal = phi2fugPhaseEq_FH_LM(params,phiFinal);
    else 
        fugFinal = 0;
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

end

