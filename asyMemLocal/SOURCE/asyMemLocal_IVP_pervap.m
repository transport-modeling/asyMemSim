%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA_IVP(z,stateVar,diffVar,localCompFlux,params) %
% Description: RHS function for solveing of initial value problem         %
% Input:       z             - spacial variable to integrate over (um)    %
%              stateVar      - (ODE, FH) n+1 dimensional vector of volume %
%                                fractions in membrane phase              %
%                              (DAE, FH-LM or DSM) 2*n+1 dimensinal vec   %
%                                of n+1 volume fractions and n fugacities %
%                                of membrane phase                        %
%              diffVar       - (ODE) diffVar = 0 (n/a)                    %
%                              (DAE) 2*n+1  dimensinal derivative vcector %
%                                of n+1 volume fractions and n fugacities %
%                                of membrane phase                        %
%              localCompFlux - n+1 dimensional vector of support layer    %
%                                compositions and total local mem flux    %
%              params        - struct of system parameters                %
%                                     (see dataBank function for specs)   %
% Output:      funIVP        - function value for integrator solver       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
function [funIVP] = asyMemLocal_IVP_pervap(z,stateVar,diffVar,localCompFlux,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters and modify as needed
    diffs = params.diffs;
    if params.memPhaseModel == 1 || params.memPhaseModel == 3
        chis = params.chis; %note chi_ji = chi_ij
    elseif params.memPhaseModel == 2
        chis = zeros(params.n+1);
    end
    n = params.n;
    HanSolParam = params.HanSolParam;
    T = params.T;
    R = params.R;
    psat = params.psat;
    Vs = params.Vs;
    Pu = params.Pu; 
    funIVP = zeros(n+1,1);
    [chis,diffs] = correlationEval(stateVar,diffs,chis,params);
    params.diffs = diffs;
    params.chis = chis;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate IVP residual functions
    if params.memPhaseModel == 1
      %solving for volume fluxes and mol fractions
        [B,invGam] = matrixEvalGammaB_FH(stateVar,params);
        if params.memPhaseModel_FicksOG == 1
            funIVP(1:n)=-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
            funIVP(n+1)=-sum(funIVP(1:n)); 
        else
            if params.noThermoCoupling == 1
                funIVP(1:n)=-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1)=-sum(funIVP(1:n));
            else
                funIVP(1:n)=-invGam*B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1)=-sum(funIVP(1:n));
            end
        end
    elseif params.memPhaseModel == 2
      %solving for volume fluxes and mol fractions
        [B,invGam] = matrixEvalGammaB_DSM(stateVar,params);
        if params.memPhaseModel_FicksOG == 1
            funIVP(1:n) = -diffVar(1:n)-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
            funIVP(n+1) = sum(diffVar(1:n+1));
            funIVP(n+2:n+1+n) = DAEevalDSM_RHS(stateVar,params); %DSM
        else
            if params.noThermoCoupling == 1
                funIVP(1:n) = -diffVar(1:n)-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1) = sum(diffVar(1:n+1));
                funIVP(n+2:n+1+n) = DAEevalDSM_RHS(stateVar,params); %DSM
            else
                funIVP(1:n) = -diffVar(1:n)-invGam*B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1) = sum(diffVar(1:n+1));
                funIVP(n+2:n+1+n) = DAEevalDSM_RHS(stateVar,params); %DSM
            end
        end
    elseif params.memPhaseModel == 3
      %solving for volume fluxes and mol fractions
        [B,invGam] = matrixEvalGammaB_FH_LM(stateVar,params);
        if params.memPhaseModel_FicksOG == 1
            funIVP(1:n) = -diffVar(1:n)-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
            funIVP(n+1) = sum(diffVar(1:n+1));
            funIVP(n+2:n+1+n) = DAEevalFH_LM_RHS(stateVar,params); %FH-DSM    
        else
            if params.noThermoCoupling == 1
                funIVP(1:n) = -diffVar(1:n)-B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1) = sum(diffVar(1:n+1));
                funIVP(n+2:n+1+n) = DAEevalFH_LM_RHS(stateVar,params); %FH-DSM 
            else
                funIVP(1:n) = -diffVar(1:n)-invGam*B*localCompFlux(1:n)*1000/3600; %note correct convertsion from LMH to UM3/UM2/s
                funIVP(n+1) = sum(diffVar(1:n+1));
                funIVP(n+2:n+1+n) = DAEevalFH_LM_RHS(stateVar,params); %FH-DSM 
            end
        end
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

end


    
   
    
    
     
            

    

            
    
            
        
    
        


