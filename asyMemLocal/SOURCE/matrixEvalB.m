%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    matrixEvalB(stateVar,params)                              %
% Description: Evaluate B matrix for MS model                             %
% Input:       stateVar - 2*n+1 dimensinal vec of n+1 volume fractions    %
%                           and n fugacities of membrane phase            %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
% Output:      B        - n x n matrix for MS eqns (diffusional terms)    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B] = matrixEvalB(stateVar,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    diffs = params.diffs;
    n = params.n;
    B = zeros(n);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies
%     for p = 1:n
%         for q = 1:n
%             if p==q
%                 if params.diffModel == 1
%                     B(p,q)=stateVar(n+1)/diffs(p,n+1);
%                 else
%                     B(p,q)=(sum(stateVar(1:n)./(diffs(p,1:n).')))-stateVar(p)/diffs(p,p)+stateVar(n+1)/diffs(p,n+1);   
%                 end   
%             else  %...p!=q
%                 if params.diffModel == 1
%                     B(p,q)=0;
%                 elseif params.diffModel == 4 && diffs(p,q)>=1E5
%                     B(p,q)=0;
%                 else
%                     B(p,q)=-stateVar(p)/diffs(p,q);
%                 end
%             end
%         end
%     end
    B = diag(stateVar(n+1)./diffs(:,n+1)+(1./diffs(:,1:n))*stateVar(1:n))-diag(stateVar(1:n))*(1./diffs(:,1:n));
%------------------------------------------------------------------------------------------------------------------------------------% 


%% tests
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies
% options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000);
% stateVar(n+2:end) = fsolve(@(fss)DAEevalFH_LM_RHS([stateVar(1:n+1);fss],params),ones(n,1)*0.002,options);
% fs = stateVar(n+2:end);
% phiFH(1:n,1) = (stateVar(1:n)-stateVar(n+1).*Ch.*bs.*fs/(1+sum(bs.*fs)));
% phiFH(n+1,1) = 1-sum(phiFH(1:n));
%     for p = 1:n
%         for q = 1:n
%             if p==q
%                 if params.diffModel == 1
%                     B(p,q)=stateVar(n+1)/diffs(p,n+1);
%                 else
%                     B(p,q)=(sum(stateVar(1:n)./(diffs(p,1:n).')))-stateVar(p)/diffs(p,p)+stateVar(n+1)/diffs(p,n+1);   
%                 end   
%                 dgdphiFH(p,q)=(1-phiFH(p)*(1+Vs(p)/Vs(q)*(phiFH(1:p-1).'*(chis(1:p-1,p)))))/phiFH(p)...
%                     +Vs(p)*(chis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-chis(n+1,p)*phiFH(p)/Vs(p))+Vs(p)/Vs(n+1)-chis(p,1+p:n)*phiFH(1+p:n)...
%                     -Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))-chis(p,n+1)*(sum(phiFH(1:n))-phiFH(p))-2*chis(p,n+1)*phiFH(n+1);
%                 dhdf(p,q) = -((1+sum(bs.*fs))*Ch(p)*bs(p)-Ch(p)*bs(p)*fs(p)*bs(p))/((1+sum(bs.*fs))^2)*stateVar(n+1);
%             else  %...p!=q
%                 if params.diffModel == 1
%                     B(p,q)=0;
%                 elseif params.diffModel == 4 && diffs(p,q)>=1E5
%                     B(p,q)=0;
%                 else
%                     B(p,q)=-stateVar(p)/diffs(p,q);
%                 end
%                 if p<q
%                     dgdphiFH(p,q)=(-Vs(p)/Vs(q)+Vs(p)*(chis(1:p-1,p).'*(phiFH(1:p-1)./Vs(1:p-1)))...
%                     +chis(p,p+1:n)*phiFH(p+1:n)-chis(p,q)*phiFH(q)+chis(p,q)*(sum(phiFH(1:n))-phiFH(p)...
%                     -phiFH(q))+2*phiFH(q)*chis(p,q)-Vs(p)*(chis(1:q-1,q).'*(phiFH(1:q-1)./Vs(1:q-1))-phiFH(p)*...
%                     chis(p,q)/Vs(p))-Vs(p)*(chis(q,q+1:n)*phiFH(q+1:n))/Vs(q))...
%                     +Vs(p)/Vs(n+1)-chis(p,1+p:n)*phiFH(1+p:n)...
%                     +chis(p,q)*phiFH(q)-Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))-chis(p,n+1)*(sum(phiFH(1:n))-phiFH(p)-phiFH(q))...
%                     +Vs(p)*(chis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-chis(n+1,p)*phiFH(p)/Vs(p)-chis(n+1,q)*phiFH(q)/Vs(q))+...
%                     (chis(p,n+1)+chis(p,q)-Vs(p)/Vs(q)*chis(q,n+1))*(phiFH(n+1)-phiFH(q))-2*phiFH(n+1)*chis(p,n+1);
%                 else %p>q
%                     dgdphiFH(p,q)=(-Vs(p)/Vs(q)+Vs(p)*(chis(1:p-1,p).'*(phiFH(1:p-1)./Vs(1:p-1))-chis(q,p)*...
%                     phiFH(q)/Vs(q))+chis(p,p+1:n)*phiFH(p+1:n)+chis(q,p)*Vs(p)*(sum(phiFH(1:n))-phiFH(p)...
%                     -phiFH(q))/Vs(q)+2*phiFH(q)*chis(p,q)*Vs(p)/Vs(q)-Vs(p)*(chis(1:q-1,q).'*(phiFH(1:q-1)./Vs(1:q-1))...
%                     )-Vs(p)*(chis(q,q+1:n)*phiFH(q+1:n))/Vs(q))...
%                     +Vs(p)/Vs(n+1)-chis(p,1+p:n)*phiFH(1+p:n)...
%                     -Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(chis(1:p-1,p))+Vs(p)*phiFH(q)*chis(q,p)/Vs(q)-chis(p,n+1)*(sum(phiFH(1:n))-phiFH(p)-phiFH(q))...
%                     +Vs(p)*(chis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-chis(n+1,p)*phiFH(p)/Vs(p)-chis(n+1,q)*phiFH(q)/Vs(q))...
%                     +(chis(n+1,p)+chis(p,q)*Vs(p)/Vs(q)-Vs(p)/Vs(q)*chis(q,n+1))*(phiFH(n+1)-phiFH(q))-2*phiFH(n+1)*chis(p,n+1);
%                 end
%                 dhdf(p,q) = Ch(p)*bs(p)*fs(p)*(bs(q)/((1+sum(bs.*fs))^2))*stateVar(n+1);
%             end
%         end
%     end
% 
% %verify math by comparing against FD approx    
% invdfdphi_OG = (-inv(eye(n)+ones(n).*Ch.*bs.*fs/(1+sum(bs.*fs)))*...
%         (dhdf-inv(dgdphiFH)*diag(1./fs)));
% %     invdfdphi = (-inv(eye(n)+ones(n).*Ch.*bs.*fs/(1+sum(bs.*fs)))*...
% %         (dhdf-inv(dgdphiFH)*diag(1./fs)));
%     for i = 1:n
%         stateVarim(1:n+1,1) = stateVar(1:n+1,1)+[zeros(i-1,1);-1E-5;zeros(n-(i),1);1E-5];
%         options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000);
%         fim = fsolve(@(fss)DAEevalFH_LM_RHS([stateVarim(1:n+1);fss],params),ones(n,1)*0.002,options);
%         fios = stateVar(n+2:end);
%         dfdphitest(:,i) = (fios-fim)./1E-5;    
%     end
%     invdfdphitest = inv(dfdphitest);
%     0
end

    
   
   