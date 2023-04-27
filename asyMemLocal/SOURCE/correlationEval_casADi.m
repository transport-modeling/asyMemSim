%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    correlationEval(stateVar,diffs,chis,params)                %
% Description: Evaluate diffusional and sorption correlations based on    %
%                stateVar values and component properties.                %
% Input:       stateVar    - (ODE, FH) n+1 dimensional vector of volume   %
%                               fractions in membrane phase               %
%                            (DAE, FH-LM or DSM) 2*n+1 dimensinal vec     %
%                              of n+1 volume fractions and n fugacities   %
%                              of membrane phase                          %
%              diffs       - n x n+1 dimensinal matrix of volume based MS %
%                              diffusivities (um^2/s)                     %
%              chis        - n+1 x n+1 dimensional matrix of FH or FH-LM  %
%                              chi parameters                             %
%              params      - struct of system parameters                  %
%                              (see dataBank function for specs)          %
%              solverState - n+1 dim variable shooting point values       %
%                              and usual localCompFlux values for current % 
%                              iteration                                  %
% Output:      diffs, chis - evaluated matricies as defined above         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chis,diffs] = correlationEval_casADi(stateVar,diffs,chis,params,solverState)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    Vs = params.Vs;
    R = params.R;
    T = params.T;
    HanSolParam = params.HanSolParam;
    yf = params.yf;
    n = params.n;
    pDiffs = params.diffs;
    unitActPhis = params.unitActPhis;
    B = params.Bffv;
%     molFracStateVar = vol2molFrac_casADi(stateVar,params);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%diffusivity modify for swelling model
    if params.swlDiffModel == 2
        for i = 1:n
              diffs(i,n+1) = pDiffs(i,n+1)*exp(B(i)*(1/unitActPhis(i)-1/abs((1-stateVar(n+1)))));
        end 
    elseif params.swlDiffModel == 3
        for i = 1:n
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*params.yf);
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*params.yf.*Vs(1:n)/sum(Vs(1:n).*params.yf));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*molFracStateVar);
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*stateVar(1:n))/sum(stateVar(1:n));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*stateVar(n+2:end)./(molFracStateVar.*purFug))/sum(stateVar(n+2:end)./(molFracStateVar.*purFug));
            diffs(i,n+1) = sum(pDiffs(1:n,n+1))/n;
%             diffs(i,n+1) = median(pDiffs(1:n,n+1));
        end 
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%diffusivity modify for cross-term Diffusivity:
    if params.diffModel == 2
        for i = 1:n
            diffs(i,i) = 1;
            for j = i+1:n
                diffs(i,j) = Vs(j)*(diffs(j,n+1)/Vs(i))^(abs(stateVar(j))/(abs(stateVar(i))+abs(stateVar(j))))*(diffs(i,n+1)/Vs(j))^(abs(stateVar(i))/(abs(stateVar(i))+abs(stateVar(j))));
                diffs(j,i) = diffs(i,j)*Vs(i)/(Vs(j));
            end
        end
        if params.crossDiffFudge == 1
            for k = 1:length(params.crossDiffSpecs)/2
                p = params.crossDiffSpecs(2*k);
                w = params.crossDiffSpecs(2*k-1);
                diffs(w,p) = params.crossDiffVals(k);
                diffs(p,w) = diffs(w,p)*params.Vs(p)/params.Vs(w);
            end
        end
    elseif params.diffModel == 3
        for i = 1:n
            diffs(i,i) = 1;
            for j=i+1:n
                diffs(i,j) = Vs(j)*((diffs(j,n+1)/Vs(i))*(abs(stateVar(j))/(abs(stateVar(i))+abs(stateVar(j))))+(diffs(i,n+1)/Vs(j))*(abs(stateVar(i))/(abs(stateVar(i))+abs(stateVar(j)))));
                diffs(j,i) = diffs(i,j)*Vs(i)/(Vs(j));
            end
        end
        if params.crossDiffFudge == 1
            for k = 1:length(params.crossDiffSpecs)/2
                p = params.crossDiffSpecs(2*k);
                w = params.crossDiffSpecs(2*k-1);
                diffs(w,p) = params.crossDiffVals(k);
                diffs(p,w) = diffs(w,p)*params.Vs(p)/params.Vs(w);
            end
        end
    elseif params.diffModel == 4
        diffs(:,1:n) = pDiffs(:,1:n);
        
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------%     
%cross-chi-parm calc
    for m = 1:n
        for k = m+1:n
            if params.pervapMode == 0
                chis(m,k) = (sum(Vs(1:n).*(yf+solverState(1:n))./2))/(R*T/9.86)*((HanSolParam(m,1)-HanSolParam(k,1))^2+...
                    0.25*(HanSolParam(m,2)-HanSolParam(k,2))^2+0.25*(HanSolParam(m,3)-HanSolParam(k,3))^2);
            elseif params.pervapMode == 1
                molFrac = (solverState(1:n)/sum(solverState(1:n)));
                chis(m,k) = (sum(Vs(1:n).*(yf+molFrac)./2))/(R*T/9.86)*((HanSolParam(m,1)-HanSolParam(k,1))^2+...
                    0.25*(HanSolParam(m,2)-HanSolParam(k,2))^2+0.25*(HanSolParam(m,3)-HanSolParam(k,3))^2);
            end
            chis(k,m) = chis(m,k);
        end
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

end

