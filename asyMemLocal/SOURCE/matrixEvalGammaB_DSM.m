%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    matrixEvalGammaB_DSM(stateVar,params)                      %
% Description: Evaluate B and invGamma matricies for MS model using dual  %
%                -sorption model (DSM).                                   %
% Input:       stateVar - 2*n+1 dimensinal vec of n+1 volume fractions    %
%                           and n fugacities of membrane phase            %
%              params   - struct of system parameters                     %
%                           (see dataBank function for specs)             %
% Output:      B        - n x n matrix for MS eqns (diffusional terms)    %
%              invGam   - n x n matrix for MS eqns (thermodynamic terms)  % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [B,invGam] = matrixEvalGammaB_DSM(stateVar,params) 

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    n = params.n;
    phis = stateVar(1:n+1);
    ks = params.ks;
    bs = params.bs;
    Ch = params.Ch;
    diffs = params.diffs;
    fs = stateVar(n+2:end);
    B = zeros(n);
%    dgdphi_inv = eye(n);
    dgdf = zeros(n);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies
    for p = 1:n
        for q = 1:n
            if p==q
                if params.diffModel == 1
                    B(p,q) = phis(n+1)/diffs(p,n+1);
                else
                    B(p,q) = (sum(phis(1:n)./(diffs(p,1:n).')))-phis(p)/diffs(p,p)+phis(n+1)/diffs(p,n+1);   
                end   
                dgdf(p,q) = -(ks(p)+((1+sum(bs.*fs))*Ch(p)*bs(p)-Ch(p)*bs(p)*fs(p)*bs(p))/((1+sum(bs.*fs))^2))*stateVar(n+1);
            else  %...p!=q
                if params.diffModel == 1
                    B(p,q) = 0;
                elseif params.diffModel == 4 && diffs(p,q)>=1E5
                    B(p,q) = 0;
                else
                    B(p,q) = -phis(p)/diffs(p,q);
                end
                dgdf(p,q) = Ch(p)*bs(p)*fs(p)*(bs(q)/((1+sum(bs.*fs))^2))*stateVar(n+1);
            end
        end
    end
    invGam = -inv(eye(n)+ones(n).*(ks.*fs + Ch.*bs.*fs/(1+sum(bs.*fs))))*dgdf*diag(fs./phis(1:n));
%------------------------------------------------------------------------------------------------------------------------------------% 
%% tests

%verify math by comparing against FD approx   
% options = optimoptions(@fsolve,'Display','off','MaxFunctionEvaluations',5000,'MaxIterations',1000);
% stateVar(n+2:end) = fsolve(@(fss)DAEevalDSM_RHS([stateVar(1:n+1);fss],params),ones(n,1)*0.2,options);
% fs = stateVar(n+2:end);
%     for p = 1:n
%         for q = 1:n
%             if p==q
%                 if params.diffModel == 1
%                     B(p,q) = phis(n+1)/diffs(p,n+1);
%                 else
%                     B(p,q) = (sum(phis(1:n)./(diffs(p,1:n).')))-phis(p)/diffs(p,p)+phis(n+1)/diffs(p,n+1);   
%                 end   
%                 dgdf(p,q) = -(ks(p)+((1+sum(bs.*fs))*Ch(p)*bs(p)-Ch(p)*bs(p)*fs(p)*bs(p))/((1+sum(bs.*fs))^2))*stateVar(n+1);
%             else  %...p!=q
%                 if params.diffModel == 1
%                     B(p,q) = 0;
%                 elseif params.diffModel == 4 && diffs(p,q)>=1E5
%                     B(p,q) = 0;
%                 else
%                     B(p,q) = -phis(p)/diffs(p,q);
%                 end
%                 dgdf(p,q) = Ch(p)*bs(p)*fs(p)*(bs(q)/((1+sum(bs.*fs))^2))*stateVar(n+1);
%             end
%         end
%     end
% invdfdphi = -inv(eye(n)+ones(n).*(ks.*fs + Ch.*bs.*fs/(1+sum(bs.*fs))))*dgdf;
%     for i = 1:n
%         stateVarim(1:n+1,1) = stateVar(1:n+1,1)+[zeros(i-1,1);-1E-6;zeros(n-(i),1);1E-6];
%         options = optimoptions(@fsolve,'Display','iter','MaxFunctionEvaluations',5000,'MaxIterations',1000);
%         fim = fsolve(@(fss)DAEevalDSM_RHS([stateVarim(1:n+1);fss],params),ones(n,1)*0.2,options);
%         fios = stateVar(n+2:end);
%         dfdphitest(:,i) = (fios-fim)./1E-6;  
% 
%     end
%     invdfdphitest = inv(dfdphitest);
%     0;
    
end
    
    

    
    

    
    
    