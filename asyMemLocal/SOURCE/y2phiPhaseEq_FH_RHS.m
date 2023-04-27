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

function [funRHS] = y2phiPhaseEq_FH_RHS(phis,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    chis = params.chis;
    diffs = params.diffs;
    ys = params.yf;
    n = params.n;
    Vs = params.Vs;
    Pu = params.Pu;
    psat = params.psat;
    R = params.R;
    T = params.T;
    W_phiFH = zeros(n,1);
    [chis,~] = correlationEval(phis,diffs,chis,params);
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
%evaluate RHS of FH
%     for i = 1:n
%         W_phiFH(i) = (1-phis(i))-Vs(i)*(sum(phis./Vs)-phis(i)/Vs(i))+Vs(i)*(sum(phis(1:i-1)...
%             .*(chis(1:i-1,i))./Vs(1:i-1)))*(sum(phis)-phis(i))+sum(phis(1+i:n+1).*((chis(i,1+i:n+1)).'))*(sum(phis)-phis(i));        
%         for j = 1:n
%                 if j == i, continue, end
%             W_phiFH(i) = W_phiFH(i)-Vs(i)/Vs(j)*sum(((chis(j,j+1:n+1)).').*phis(j+1:n+1))*phis(j);
%         end  
%     end
    oneMinus = ones(n,1)-phis(1:n);
    W_phiFH = oneMinus+(FH_VFV*(-eye(n+1)-diag([upTriChis*phis;0]))...
        +diag(oneMinus)*((FH_VLV+upTriChis)))*phis;
    funRHS = ys-phis(1:n).*exp(W_phiFH);
    funRHS(n+1) = 1-sum(phis);
%------------------------------------------------------------------------------------------------------------------------------------% 

end
   
    