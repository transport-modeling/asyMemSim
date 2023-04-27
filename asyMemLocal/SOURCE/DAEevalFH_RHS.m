%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    evalFH_RHS(stateVar,params)                                %
% Description: Evaluate RHS of Flory-Huggings equation.                   %
% Input:       stateVar - (ODE, FH) n+1 dimensional vector of volume      %
%                           fractions in membrane phase                   %
%                         (DAE, FH-LM or DSM) 2*n+1 dimensinal vector     %
%                           of n+1 volume fractions and n fugacities      %
%                           of membrane phase                             %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
% Output:      FH_RHS   - nonlinear function value of FH RHS equation     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [FH_RHS] = DAEevalFH_RHS(stateVar,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters, intialize, and solve variables for as needed
    chis = params.chis;
    Vs = params.Vs;
    n = params.n;   
    fs = stateVar(n+2:end);
    Pu = params.Pu;
    psat = params.psat;
    R = params.R;
    T = params.T;
    phiFH(1:n+1,1) = stateVar(1:n+1);
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
    oneMinus = ones(n,1)-phiFH(1:n);
    W_phiFH = oneMinus+(FH_VFV*(-eye(n+1)-diag([upTriChis*phiFH;0]))...
        +diag(oneMinus)*((FH_VLV+upTriChis)))*phiFH;
%     W_phiFH = oneMinus+([-1-upTriChis*phiFH;-1].'.*FH_VFV+oneMinus.*(FH_VLV+upTriChis))*phiFH;   %need to testes
%      FH_RHS = phiFH(1:n).*exp(W_phiFH)-fs./(psat.*exp(Vs(1:n).*(Pu-psat*0.00131)/(R*T)));
   FH_RHS = phiFH(1:n).*exp(W_phiFH)-fs; %if assuming fs = activity
%------------------------------------------------------------------------------------------------------------------------------------% 

end

