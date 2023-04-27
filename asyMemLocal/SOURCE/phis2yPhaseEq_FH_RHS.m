%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    phi2yPhaseEq_FH_RHS(phis,params)                           %
% Description: Evaluate RHS function for solving of feed side membrane    %
%              phase vol frac based on Flory-Huggins (FH) sorption model. % 
% Input:       params - struct of system parameters                       %
%                         (see dataBank function for specs)               %
%              phis   - n+1 dimensional vector of membrane phase volume   %
%                         fractions                                       %
% Output:      ys     - mole fractions in equalibrium with input phis     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [ys] = phis2yPhaseEq_FH_RHS(phis,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    Vs = params.Vs;
    chis = params.chis;
    psat = params.psat;
    Pu = params.Pu;
    Pd = params.Pd;
    R = params.R;
    T = params.T;
    n = params.n;  
    W_phiFH = zeros(n,1);
    upTriChis = params.upTriChis;
    FH_VFV = params.FH_VFV;
    FH_VLV = params.FH_VLV;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate FH RHS
%     for i = 1:n
%         W_phiFH(i) = (1-phiFH(i))...     
%         +chis(i,n+1)*phiFH(n+1)*phiFH(n+1);
%             -Vs(i)*(sum(phiFH./Vs)-phiFH(i)/Vs(i))+Vs(i)*(sum(phiFH(1:i-1)...
%             .*(chis(1:i-1,i))./Vs(1:i-1)))*(sum(phiFH)-phiFH(i))+sum(phiFH(1+i:n+1).*((chis(i,1+i:n+1)).'))*(sum(phiFH)-phiFH(i));   
%         for j = 1:n
%             if j == i, continue, end
%             W_phiFH(i) = W_phiFH(i)-Vs(i)/Vs(j)*sum(((chis(j,j+1:n+1)).').*phiFH(j+1:n+1))*phiFH(j);
%         end  
%     end
    oneMinus = ones(n,1)-phis(1:n);
    W_phiFH = oneMinus+(FH_VFV*(-eye(n+1)-diag([upTriChis*phis;0]))...
        +diag(oneMinus)*((FH_VLV+upTriChis)))*phis;
  %function output
  ys = phis(1:n).*exp(W_phiFH)./exp(-Vs(1:n).*(Pu-Pd)/(R*T));

%------------------------------------------------------------------------------------------------------------------------------------% 

end
