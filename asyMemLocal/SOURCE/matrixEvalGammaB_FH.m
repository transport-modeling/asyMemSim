%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    matrixEvalGammaB_FH(stateVar,params)                       %
% Description: Evaluate B and invGamma matricies for MS model using Flory %
%                -Huggins sorption model (FH).                            %
% Input:       stateVar - n+1 dimensional vector of volume fractions in   %
%                           membrane phase                                %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
% Output:      B        - n x n matrix for MS eqns (diffusional terms)    %
%              invGam   - n x n matrix for MS eqns (thermodynamic terms)  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,invLam]= matrixEvalGammaB_FH(stateVar,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    diffs = params.diffs;
    n = params.n;
    Vs = params.Vs;
    chis = params.chis;
    [chis,diffs] = correlationEval(stateVar,diffs,chis,params);
    B = zeros(n);
    gamma = zeros(n);
    upTriChis = params.upTriChis;
    FH_VFV = params.FH_VFV;
    FH_VLV = params.FH_VLV;
%     stateVar = stateVar./1000;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies
%     for p = 1:n
%         for q = 1:n
%             if p==q
%                if params.diffModel == 1
%                     B(p,q)=stateVar(n+1)/diffs(p,n+1);
%                 else
%                     B(p,q)=(sum(stateVar(1:n)./(diffs(p,1:n).')))-stateVar(p)/diffs(p,p)+stateVar(n+1)/diffs(p,n+1);   
%                 end   
%                 gamma(p,q)= stateVar(p)*((1-stateVar(p)*(1+Vs(p)/Vs(q)*(stateVar(1:p-1).'*(chis(1:p-1,p)))))/stateVar(p)...
%                     +Vs(p)*(chis(n+1,1:n)*(stateVar(1:n)./Vs(1:n))-chis(n+1,p)*stateVar(p)/Vs(p))+Vs(p)/Vs(n+1)-chis(p,1+p:n)*stateVar(1+p:n)...
%                     -Vs(p)*(stateVar(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))-chis(p,n+1)*(sum(stateVar(1:n))-stateVar(p))-2*chis(p,n+1)*stateVar(n+1));
%             else  %...p!=q
%                 if params.diffModel == 1
%                     B(p,q)=0;
%                 elseif params.diffModel == 4 && diffs(p,q)>=1E5
%                     B(p,q)=0;
%                 else
%                     B(p,q)=-stateVar(p)/diffs(p,q);
%                 end
%                 if p<q
%                     gamma(p,q)=stateVar(p)*((-Vs(p)/Vs(q)+Vs(p)*(chis(1:p-1,p).'*(stateVar(1:p-1)./Vs(1:p-1)))...
%                     +chis(p,p+1:n)*stateVar(p+1:n)-chis(p,q)*stateVar(q)+chis(p,q)*(sum(stateVar(1:n))-stateVar(p)...
%                     -stateVar(q))+2*stateVar(q)*chis(p,q)-Vs(p)*(chis(1:q-1,q).'*(stateVar(1:q-1)./Vs(1:q-1))-stateVar(p)*...
%                     chis(p,q)/Vs(p))-Vs(p)*(chis(q,q+1:n)*stateVar(q+1:n))/Vs(q))...
%                     +Vs(p)/Vs(n+1)-chis(p,1+p:n)*stateVar(1+p:n)...
%                     +chis(p,q)*stateVar(q)-Vs(p)*(stateVar(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))-chis(p,n+1)*(sum(stateVar(1:n))-stateVar(p)-stateVar(q))...
%                     +Vs(p)*(chis(n+1,1:n)*(stateVar(1:n)./Vs(1:n))-chis(n+1,p)*stateVar(p)/Vs(p)-chis(n+1,q)*stateVar(q)/Vs(q))+...
%                     (chis(p,n+1)+chis(p,q)-Vs(p)/Vs(q)*chis(q,n+1))*(stateVar(n+1)-stateVar(q))-2*stateVar(n+1)*chis(p,n+1));
%                 else %p>q
%                     gamma(p,q)=stateVar(p)*((-Vs(p)/Vs(q)+Vs(p)*(chis(1:p-1,p).'*(stateVar(1:p-1)./Vs(1:p-1))-chis(q,p)*...
%                     stateVar(q)/Vs(q))+chis(p,p+1:n)*stateVar(p+1:n)+chis(q,p)*Vs(p)*(sum(stateVar(1:n))-stateVar(p)...
%                     -stateVar(q))/Vs(q)+2*stateVar(q)*chis(p,q)*Vs(p)/Vs(q)-Vs(p)*(chis(1:q-1,q).'*(stateVar(1:q-1)./Vs(1:q-1))...
%                     )-Vs(p)*(chis(q,q+1:n)*stateVar(q+1:n))/Vs(q))...
%                     +Vs(p)/Vs(n+1)-chis(p,1+p:n)*stateVar(1+p:n)...
%                     -Vs(p)*(stateVar(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))+Vs(p)*stateVar(q)*chis(q,p)/Vs(q)-chis(p,n+1)*(sum(stateVar(1:n))-stateVar(p)-stateVar(q))...
%                     +Vs(p)*(chis(n+1,1:n)*(stateVar(1:n)./Vs(1:n))-chis(n+1,p)*stateVar(p)/Vs(p)-chis(n+1,q)*stateVar(q)/Vs(q))...
%                     +(chis(n+1,p)+chis(p,q)*Vs(p)/Vs(q)-Vs(p)/Vs(q)*chis(q,n+1))*(stateVar(n+1)-stateVar(q))-2*stateVar(n+1)*chis(p,n+1));
%                 end
%             end
%         end
%     end
    B = diag(stateVar(n+1)./diffs(:,n+1)+(1./diffs(:,1:n))*stateVar(1:n))-diag(stateVar(1:n))*(1./diffs(:,1:n));
    oneMinus = diag(ones(n,1)-stateVar(1:n));
    eyeBones = [eye(n);-ones(1,n)];
    phisMatr = stateVar(1:n+1).*ones(n+1,n);
    gamma = stateVar(1:n).*(diag(1./stateVar(1:n))-eye(n)-FH_VFV*(eyeBones+diag(-chis(:,n+1))*phisMatr+diag([upTriChis*stateVar(1:n+1);0])*eyeBones)...
        -diag(diag((FH_VLV+upTriChis)*phisMatr))+oneMinus*(FH_VLV+upTriChis)*eyeBones);
%     gamma_Vec_1 = stateVar(1:n).*([1./stateVar(1);0;0;0;0;0;0;0;0]-[1;0;0;0;0;0;0;0;0]-FH_VFV*([1;0;0;0;0;0;0;0;0;-1]+diag(-chis(:,n+1))*phisMatr(:,1)+diag([upTriChis*stateVar(1:n+1);0])*[1;0;0;0;0;0;0;0;0;-1])...
%         -[[1;0;0;0;0;0;0;0;0],zeros(9,8)]*(FH_VLV+upTriChis)*phisMatr(:,1)+oneMinus*(FH_VLV+upTriChis)*[1;0;0;0;0;0;0;0;0;-1]);
%     gamma_Vec_2 = stateVar(1:n).*([0;1./stateVar(2);0;0;0;0;0;0;0]-[0;1;0;0;0;0;0;0;0]-FH_VFV*([0;1;0;0;0;0;0;0;0;-1]+diag(-chis(:,n+1))*phisMatr(:,2)+diag([upTriChis*stateVar(1:n+1);0])*[0;1;0;0;0;0;0;0;0;-1])...
%         -[[0;0;0;0;0;0;0;0;0],[0;1;0;0;0;0;0;0;0],zeros(9,7)]*(FH_VLV+upTriChis)*phisMatr(:,2)+oneMinus*(FH_VLV+upTriChis)*[0;1;0;0;0;0;0;0;0;-1]);
%     gamma_Vec_3 = stateVar(1:n).*([0;0;1./stateVar(3);0;0;0;0;0;0]-[0;0;1;0;0;0;0;0;0]-FH_VFV*([0;0;1;0;0;0;0;0;0;-1]+diag(-chis(:,n+1))*phisMatr(:,3)+diag([upTriChis*stateVar(1:n+1);0])*[0;0;1;0;0;0;0;0;0;-1])...
%         -[[0;0;0;0;0;0;0;0;0],[0;0;0;0;0;0;0;0;0],[0;0;1;0;0;0;0;0;0],zeros(9,6)]*(FH_VLV+upTriChis)*phisMatr(:,3)+oneMinus*(FH_VLV+upTriChis)*[0;0;1;0;0;0;0;0;0;-1]);
    invLam = inv(gamma);  %FH type sopriton throughout mem
   
  
%------------------------------------------------------------------------------------------------------------------------------------% 
%   for i = 1:n
%         stateVarim(1:n+1,1) = stateVar(1:n+1,1)+[zeros(i-1,1);-1E-5;zeros(n-(i),1);1E-5];
%         options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000,'FunctionTolerance',10E-8);
%         fim = fsolve(@(fss)DAEevalFH_RHS([stateVarim(1:n+1);fss],params),ones(n,1)*0.002,options);
%         fios = fsolve(@(fss)DAEevalFH_RHS([stateVar(1:n+1);fss],params),ones(n,1)*0.002,options);
%         gammaTest(:,i) = stateVar(1:n).*(log(fios)-log(fim))./1E-5;
%     end
%     A = (gammaTest-gamma);
%     B = (gamma_Vec-gamma);
%     C = (gamma_Vec-gammaTest);
end


    
    
    