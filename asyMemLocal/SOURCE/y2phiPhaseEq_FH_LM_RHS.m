%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    y2phiPhaseEq_FH_RHS(phis,params)                           %
% Description: Evaluate RHS function for solving of feed side membrane    %
%              phase vol frac based on Flory-Huggins (FH) sorption model. % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
%              phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
% Output:      funRHS - function values vector for FH equations           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FH_RHS] = y2phiPhaseEq_FH_LM_RHS(stateVar,params)
%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters, intialize, and solve variables for as needed
    Ch = params.Ch;
    bs = params.bs;
    chis = params.chis;
    Vs = params.Vs;
    n = params.n;   
    fs = stateVar(n+2:end);
    Pu = params.Pu;
    psat = params.psat;
    R = params.R;
    T = params.T;
    phiFH(1:n,1) = (stateVar(1:n)-stateVar(n+1).*Ch.*bs.*fs/(1+sum(bs.*fs)));
    phiFH(n+1,1) = 1-sum(phiFH(1:n));
    W_phiFH = zeros(n,1);
    [chis,~] = correlationEval(stateVar,params.diffs,chis,params);
    frameMatrix = ones(n,n+1)-[diag(ones(n,1)),zeros(n,1)];
    upTriChis = triu(chis(1:n,:),1);
    diagVsInv = diag(1./Vs);
    diagVs = diag(Vs(1:n));
    params.upTriChis = upTriChis;
    params.FH_VFV = (diagVs*(frameMatrix)*diagVsInv);
    params.FH_VLV = diagVs*tril(chis(1:n,:),-1)*diagVsInv;
    upTriChis = params.upTriChis;
    FH_VFV = params.FH_VFV;
    FH_VLV = params.FH_VLV;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate FH RHS
%     for i = 1:n
%         W_phiFH(i) = (1-phiFH(i))-Vs(i)*(sum(phiFH./Vs)-phiFH(i)/Vs(i))+Vs(i)*(sum(phiFH(1:i-1)...
%             .*(chis(1:i-1,i))./Vs(1:i-1)))*(sum(phiFH)-phiFH(i))+sum(phiFH(1+i:n+1).*((chis(i,1+i:n+1)).'))*(sum(phiFH)-phiFH(i));        
%         for j = 1:n
%             if j == i, continue, end
%             W_phiFH(i) = W_phiFH(i)-Vs(i)/Vs(j)*sum(((chis(j,j+1:n+1)).').*phiFH(j+1:n+1))*phiFH(j);
%         end  
%     end
    oneMinus = ones(n,1)-phiFH(1:n);
    W_phiFH = oneMinus+(FH_VFV*(-eye(n+1)-diag([upTriChis*phiFH;0]))...
        +diag(oneMinus)*((FH_VLV+upTriChis)))*phiFH;
    FH_RHS = phiFH(1:n).*exp(W_phiFH)-stateVar(n+2:end)./(psat.*exp(Vs(1:n).*(Pu-psat*0.00131)/(R*T)));
    FH_RHS(n+1) = 1-sum(stateVar(1:n+1));
%------------------------------------------------------------------------------------------------------------------------------------% 
end

