%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    MSAmatrixEvalGammaB_FH-LM(stateVar,params)                 %
% Description: Evaluate B and invGamma matricies for MS model using       %
%                combined Flory-Huggins and Langmuir sorp model (FH-LM).  %
% Input:       stateVar - 2*n+1 dimensinal vec of n+1 volume fractions    %
%                           and n fugacities of membrane phase            %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
%              solverState - n+1 dim variable shooting point values       %
%                              and usual localCompFlux values for current % 
%                              iteration                                  %
% Output:      B        - n x n matrix for MS eqns (diffusional terms)    %
%              invGam   - n x n matrix for MS eqns (thermodynamic terms)  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,invGam] = matrixEvalGammaB_FH_LM_casADi(stateVar,params,solverState)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    import casadi.*
    n = params.n;
    Vs = params.Vs;
    HanSolParam = params.HanSolParam;
    phiFH = SX.sym('phi',n+1);
    dgdphiFH = SX.sym('dgdphi',n,n);
    dhdf = SX.sym('dhdf',n,n);
    B = SX.sym('B',n,n);
    sxDiffs = SX.sym('sxDiffs',n,n+1);
    sxChis = SX.sym('sxChis',n+1,n+1);
    sxChis(n+1,:) = params.chis(n+1,:);
    sxChis(:,n+1) = params.chis(:,n+1);
    if n==1
        sxDiffs(1,n+1) = params.diffs(1,n+1);
    else
        sxDiffs(:,n+1) = params.diffs(1:n,n+1);
    end
    [sxChis,sxDiffs] = correlationEval_casADi(stateVar,sxDiffs,sxChis,params,solverState);
    bs = params.bs;
    Ch = params.Ch;
    fs = stateVar(n+2:end);
    phiFH(1:n,1) = (stateVar(1:n)-stateVar(n+1).*Ch.*bs.*fs/(1+sum(bs.*fs)));
    phiFH(n+1,1) = 1-sum(phiFH(1:n));
    Pu = params.Pu;
    psat = params.psat;
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies  
    for p = 1:n
        for q = 1:n
            if p==q
                if params.diffModel == 1
                    B(p,q)=stateVar(n+1)/sxDiffs(p,n+1);
                else
                    B(p,q)=(sum(stateVar(1:n)./(sxDiffs(p,1:n).')))-stateVar(p)/sxDiffs(p,p)+stateVar(n+1)/sxDiffs(p,n+1);   
                end
                dgdphiFH(p,q)=(1-phiFH(p)*(1+Vs(p)/Vs(q)*(phiFH(1:p-1).'*(sxChis(1:p-1,p)))))/phiFH(p)...
                    +Vs(p)*(sxChis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-sxChis(n+1,p)*phiFH(p)/Vs(p))+Vs(p)/Vs(n+1)-sxChis(p,1+p:n)*phiFH(1+p:n)...
                    -Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(sxChis(1:p-1,p))-sxChis(p,n+1)*(sum(phiFH(1:n))-phiFH(p))-2*sxChis(p,n+1)*phiFH(n+1);
                dhdf(p,q) = -((1+sum(bs.*fs))*Ch(p)*bs(p)-Ch(p)*bs(p)*fs(p)*bs(p))/((1+sum(bs.*fs))^2)*stateVar(n+1);
            else  %...p!=q
                if params.diffModel == 1
                    B(p,q)=0;
                elseif params.diffModel == 4 && sxDiffs(p,q)>=1E5
                    B(p,q)=0;
                else
                    B(p,q)=-stateVar(p)/sxDiffs(p,q);
                end                
                if p<q
                    dgdphiFH(p,q)=(-Vs(p)/Vs(q)+Vs(p)*(sxChis(1:p-1,p).'*(phiFH(1:p-1)./Vs(1:p-1)))...
                    +sxChis(p,p+1:n)*phiFH(p+1:n)-sxChis(p,q)*phiFH(q)+sxChis(p,q)*(sum(phiFH(1:n))-phiFH(p)...
                    -phiFH(q))+2*phiFH(q)*sxChis(p,q)-Vs(p)*(sxChis(1:q-1,q).'*(phiFH(1:q-1)./Vs(1:q-1))-phiFH(p)*...
                    sxChis(p,q)/Vs(p))-Vs(p)*(sxChis(q,q+1:n)*phiFH(q+1:n))/Vs(q))...
                    +Vs(p)/Vs(n+1)-sxChis(p,1+p:n)*phiFH(1+p:n)...
                    +sxChis(p,q)*phiFH(q)-Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(sxChis(1:p-1,p))-sxChis(p,n+1)*(sum(phiFH(1:n))-phiFH(p)-phiFH(q))...
                    +Vs(p)*(sxChis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-sxChis(n+1,p)*phiFH(p)/Vs(p)-sxChis(n+1,q)*phiFH(q)/Vs(q))+...
                    (sxChis(p,n+1)+sxChis(p,q)-Vs(p)/Vs(q)*sxChis(q,n+1))*(phiFH(n+1)-phiFH(q))-2*phiFH(n+1)*sxChis(p,n+1);
                else %p>q
                    dgdphiFH(p,q)=(-Vs(p)/Vs(q)+Vs(p)*(sxChis(1:p-1,p).'*(phiFH(1:p-1)./Vs(1:p-1))-sxChis(q,p)*...
                    phiFH(q)/Vs(q))+sxChis(p,p+1:n)*phiFH(p+1:n)+sxChis(q,p)*Vs(p)*(sum(phiFH(1:n))-phiFH(p)...
                    -phiFH(q))/Vs(q)+2*phiFH(q)*sxChis(p,q)*Vs(p)/Vs(q)-Vs(p)*(sxChis(1:q-1,q).'*(phiFH(1:q-1)./Vs(1:q-1))...
                    )-Vs(p)*(sxChis(q,q+1:n)*phiFH(q+1:n))/Vs(q))...
                    +Vs(p)/Vs(n+1)-sxChis(p,1+p:n)*phiFH(1+p:n)...
                    -Vs(p)*(phiFH(1:p-1)./Vs(1:p-1)).'*(sxChis(1:p-1,p))+Vs(p)*phiFH(q)*sxChis(q,p)/Vs(q)-sxChis(p,n+1)*(sum(phiFH(1:n))-phiFH(p)-phiFH(q))...
                    +Vs(p)*(sxChis(n+1,1:n)*(phiFH(1:n)./Vs(1:n))-sxChis(n+1,p)*phiFH(p)/Vs(p)-sxChis(n+1,q)*phiFH(q)/Vs(q))...
                    +(sxChis(n+1,p)+sxChis(p,q)*Vs(p)/Vs(q)-Vs(p)/Vs(q)*sxChis(q,n+1))*(phiFH(n+1)-phiFH(q))-2*phiFH(n+1)*sxChis(p,n+1);
                end
                dhdf(p,q) = Ch(p)*bs(p)*fs(p)*(bs(q)/((1+sum(bs.*fs))^2))*stateVar(n+1);
            end
        end
    end
    invGam = (-inv(eye(n)+ones(n).*Ch.*bs.*fs/(1+sum(bs.*fs)))*...
        (dhdf-inv(dgdphiFH)*diag(1./fs))*diag(fs./stateVar(1:n)));
%------------------------------------------------------------------------------------------------------------------------------------% 
    
end
%% test   

%verify math by comparing against FD approx    
% dfdphi = inv(-(dhdf+diag(1-Ch)*(inv(dgdphiFH))*diag(-1./fs)));
%     for i = 1:n+1
%         fios = W(n+2:end);
%         Wim(1:n+1,1) = W(1:n+1,1)+[zeros(i-1,1);-1E-3;zeros(n+1-(i),1)];
%         options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000);
%         fim = fsolve(@(fss)RHSphitofugFHDSM(fss,Wim(1:n+1),Vs,n,Chis,Pu,bs,Ch,psat),ones(n,1)*0.000002,options);
%         dfdphitest(:,i) = (fios-fim)./1E-3;    
%     end
    

    
   
   