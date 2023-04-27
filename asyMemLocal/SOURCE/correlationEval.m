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
% Output:      diffs, chis - evaluated matricies as defined above         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [chis,diffs] = correlationEval(stateVar,diffs,chis,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    yf = params.yf;
    Vs = params.Vs;
    R = params.R;
    T = params.T;
    HanSolParam = params.HanSolParam;
    n = params.n;
    pDiffs = params.diffs;
    unitActPhis = params.unitActPhis;
    B = params.Bffv;
%     molFracStateVar = vol2molFrac(stateVar,params);
    purFug = params.purFug;
    if length(stateVar) < 2*n+1
        stateVar = [stateVar;ones(n,1)];
    end
    stateVar = stateVar./1000;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%diffusivity modify for swelling model
    if params.swlDiffModel == 2 && params.diffsDone == 0
        for i = 1:n
              diffs(i,n+1) = pDiffs(i,n+1)*exp(B(i)*(1/unitActPhis(i)-1/abs((1-stateVar(n+1)))));
        end 
    elseif params.swlDiffModel == 3
        for i = 1:n
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*params.yf);
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*params.yf.*Vs(1:n)/sum(Vs(1:n).*params.yf));
%             diffs(i,n+1) = prod((pDiffs(1:n,n+1)).^(stateVar(1:n)/(1-stateVar(n+1))));
%             if params.memPhaseModel == 1 || params.memPhaseModel == 3
%                 diffs(i,n+1) = (params.chis(i,n+1))^(-1)*prod((pDiffs(1:n,n+1).*params.chis(1:n,n+1)).^(stateVar(1:n)/(1-stateVar(n+1))));
%             else
%                 diffs(i,n+1) = (params.ks(i))^(-1)*prod((pDiffs(1:n,n+1).*params.ks(1:n)).^(stateVar(1:n)/(1-stateVar(n+1))));
%             end
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*stateVar(n+2:end)./(molFracStateVar.*purFug))/sum(stateVar(n+2:end)./(molFracStateVar.*purFug));
             diffs(i,n+1) = (Vs(i))^(-1)*prod((pDiffs(1:n,n+1).*Vs(1:n)).^(stateVar(1:n)/(1-stateVar(n+1))));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*(molFracStateVar));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*stateVar(1:n).*Vs(1:n))/sum(stateVar(1:n).*Vs(1:n));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*stateVar(1:n)/sum(stateVar(1:n)));
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1))/n;
%            diffs(i,n+1) = median(pDiffs(1:n,n+1));
%             if (params.memID == 1 && n == 9) || (params.memID == 2)
%                 compID_temp = [];
%                 aromaticsID = ["TOL","MCH","DEC","MNP","TBB","TPB","PXY","OXY"];
%                 for m = 1:length(aromaticsID)
%                     if isfield(params.compID,aromaticsID(m))
%                         compID_i = getfield(params.compID,aromaticsID(m));
%                         compID_temp = [compID_temp;compID_i];
%                     end
%                 end
%                 diffs(i,n+1) = sum(pDiffs(compID_temp,n+1))/length(compID_temp);
%             elseif (params.memID == 1 && n == 5) || (params.memID == 1 && isfield(params.compID,"IOC") ...
%                     && isfield(params.compID,"ICT")) || (params.memID == 1 && isfield(params.compID,"NOC") ...
%                     && isfield(params.compID,"ICT"))
%                 compID_temp = [];
%                 HCchainID = ["ICT","IOC","NOC"];
%                 for m = 1:length(HCchainID)
%                     if isfield(params.compID,HCchainID(m))
%                         compID_i = getfield(params.compID,HCchainID(m));
%                         compID_temp = [compID_temp;compID_i];
%                     end
%                 end
%                 diffs(i,n+1) = sum(pDiffs(compID_temp,n+1))/length(compID_temp);
%             else
%                 diffs(i,n+1) = sum(pDiffs(1:n,n+1))/(n);
%             end
%             
%             diffs(i,n+1) = sum(pDiffs(1:n,n+1).*Vs(1:n)./chis(1:n,n+1))/sum(Vs(1:n)./chis(1:n,n+1));
        end 
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%diffusivity modify for cross-term Diffusivity:
    if params.diffModel == 2 && params.diffsDone == 0
%         for i = 1:n
%             diffs(i,i) = 1;
%             for j = i+1:n
%                 diffs(i,j) = Vs(j)*(diffs(j,n+1)/Vs(i))^(abs(stateVar(j))/(abs(stateVar(i))+abs(stateVar(j))))*(diffs(i,n+1)/Vs(j))^(abs(stateVar(i))/(abs(stateVar(i))+abs(stateVar(j))));
%                 diffs(j,i) = diffs(i,j)*Vs(i)/(Vs(j));
%             end
%         end
        diffsUp = ((((diffs(1:n,n+1))).^abs(stateVar(1:n))).*(((diffs(1:n,n+1)).^abs(stateVar(1:n))).').*(((1./(Vs(1:n)).')).^(abs(stateVar(1:n)))).*(((1./Vs(1:n)))...
            .^(stateVar(1:n).'))).^(1./(abs(stateVar(1:n)).'+abs(stateVar(1:n)))).*(Vs(1:n).');
        diffsDw = diffsUp.'.*(Vs(1:n).'./(Vs(1:n)));
        diffs = [triu(diffsUp,1) + tril(diffsDw,-1)+eye(n),diffs(:,n+1)];
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
if params.chisDone == 0
    for m = 1:n
        for k = m+1:n
            chis(m,k) = ((Vs(m)*Vs(k))^0.5)/(R*T/9.86)*((HanSolParam(m,1)-HanSolParam(k,1))^2+...
                0.25*(HanSolParam(m,2)-HanSolParam(k,2))^2+0.25*(HanSolParam(m,3)-HanSolParam(k,3))^2);
            chis(k,m) = chis(m,k);
        end
    end
elseif params.chisDone == 1
    0;
end
%------------------------------------------------------------------------------------------------------------------------------------% 

end

