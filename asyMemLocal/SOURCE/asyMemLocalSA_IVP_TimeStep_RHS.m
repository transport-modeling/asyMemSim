%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA_IVP_TimeStep_RHS(z,currentPhi,prevPhi,...    %
%                diffVar,localCompFlux,params,h)                          %
%                              species, x_scale_form)                     %
% Description: Evaluate RHS function for initial value problem (IVP)      %
%                using explicit or implicit Euler's method.               %
% Input:       z             -  spatial membrane thickness variable (um)  %
%              currentPhi    -  current step spacial volume fractions     %
%              previousPhi   -  previous step spacial volume fractions    %
%              localCompFlux - n+1 dimensional vector of support layer    %
%                                compositions and total local mem flux    %
%              params        - struct of system parameters                %
%                                (see dataBank function for specs)        %
%              h             - step size of spacial interval (um)         %
% Output:      funIVP        -  func value vector for nonlinear solver    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funIVP] = asyMemLocalSA_IVP_TimeStep_RHS(z,currentPhi,prevPhi,diffVar,localCompFlux,params,h)

%------------------------------------------------------------------------------------------------------------------------------------% 
%build RHS for Newton solve of each node of discrutuzed IVP

    %Backward Euler
%    funIVP = currentPhi-prevPhi-asyMemLocalSA_IVP(z,currentPhi,diffVar,localCompFlux,params).*h;
    
    %Crank-Nicolson
    funIVP = currentPhi-prevPhi-(asyMemLocalSA_IVP(z,currentPhi,diffVar,localCompFlux,params)+...
        asyMemLocalSA_IVP(z,prevPhi,diffVar,localCompFlux,params))/2.*h;
%------------------------------------------------------------------------------------------------------------------------------------% 

end

