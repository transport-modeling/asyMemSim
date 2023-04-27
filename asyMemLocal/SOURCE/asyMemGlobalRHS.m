%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA_IVP(z,stateVar,diffVar,localCompFlux,params) %
% Description: RHS function for solveing of initial value problem         %
% Input:       x        - spacial variable to integrate over (m)          %
%              stateVar - vector of n componet volume fractions and       % 
%                           state localCompFlux                           %      
%              diffVar  - vector of n derivatives of componet volume      %
%                           fractions and derivative of localCompFlux     % 
%              params   - struct of system parameters                     %
%                                     (see dataBank function for specs)   %
% Output:      resDAE   - function value for integrator DAE solver        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [resDAE] = asyMemGlobalRHS(x,stateVar,diffVar,params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    diffs = params.diffs;
    chis = params.chis;
    n = params.n;
    bs = params.bs;
    Ch = params.Ch;
    psat = params.psat;
    T = params.T;
    R = params.R;
    Vs = params.Vs;
    Pu = params.Pu;
    yf = params.yf;
    %ks = params.ks;
    Pd = params.Pd;
    %variable perm
    params.yf = stateVar(1:n)/sum(stateVar(1:n+1));
    %cons perm (note change if Cintert change)
%     params.yf = [params.yfinit(1:n)/sum([params.yfinit;0.01])];
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve for feed BC
    phiFeed = bulk_y_to_phis_FH_solve(params);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%evaluate matricies and correlations based on state
  %variable perm
    [chis,diffs] = correlationEval([(stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n-1)+phiFeed(1:n))/2;...
    1-sum((stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n-1)+phiFeed(1:n))/2)],diffs,chis,params);
  %cons perm
%     [Chis,Diffs] = CorrelationEval(phifs,Diffs,Chis,params);
    
    params.Diffs = diffs;
    params.Chis = chis;
    %note to have FH-LM need to add algebraic equations to describe fugPerm
    
  %variable permeablity
    [invLamavg,Bavg]=BLammatrixFHcreateevalij((stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n)+phiFeed)/2,params);
  %constant permeabilty (eval @ feed) note if want constant flux then comment out line 23 and 16
%     [invLamavg,Bavg]=BLammatrixFHcreateevalij(phifs,params);
    
    width = 1; %m
    resDAE = zeros(4*(n+1)+2,1);

    params.yf = stateVar(1:n)/sum(stateVar(1:n+1));
    phiFeed = bulkytophisFHsolve(params);
%------------------------------------------------------------------------------------------------------------------------------------%     
   
%------------------------------------------------------------------------------------------------------------------------------------% 
%DAE residual evaluations
  %feed species balance
    resDAE(1:n) = diffVar(1:n)*stateVar(n+2)+diffVar(n+2)*stateVar(1:n)+width*stateVar(end-n:end-1)*stateVar(end)/sum(Vs(1:n).*stateVar(end-n:end-1))*1000; 
    resDAE(n+1) = diffVar(n+1)*stateVar(n+2)+diffVar(n+2)*stateVar(n+1);
    resDAE(n+2) = diffVar(n+2)*sum(stateVar(1:n+1))+width*stateVar(end)/sum(Vs(1:n).*stateVar(end-n:end-1))*1000; %total feed balance

  %perm balance
    resDAE(n+3:n+3+n-1) = diffVar(n+3:n+3+n-1)*stateVar(2*(n+1)+2)+diffVar(2*(n+1)+2)*stateVar(n+3:n+3+n-1)-width*stateVar(end-n:end-1)*stateVar(end)/sum(Vs(1:n).*...
    stateVar(end-n:end-1))*1000; 
    resDAE(n+3+n) = diffVar(n+3+n)*stateVar(2*(n+1)+2)+diffVar(2*(n+1)+2)*stateVar(n+3+n);
    resDAE(2*(n+1)+2) = diffVar(2*(n+1)+2)*sum(stateVar(n+3:n+3+n))-width*stateVar(end)/sum(Vs(1:n).*stateVar(end-n:end-1))*1000; %total perm balance

  %MS eqns diff eqns (converts everything to ummmm)
    resDAE(2*(n+1)+2+1:2*(n+1)+2+1+n-1) = stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n-1)-phiFeed(1:n)+params.lmem*invLamavg*...
    Bavg*stateVar(end-n:end-1).*Vs(1:n).*stateVar(end)./sum(Vs(1:n).*stateVar(end-n:end-1))*1000/3600; 
    resDAE(2*(n+1)+2+1+n) = stateVar(2*(n+1)+2+1+n)-phiFeed(n+1)-sum(resDAE(2*(n+1)+2+1:2*(n+1)+2+1+n-1)-...
    (stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n-1)-phiFeed(1:n)));

  %LFM nonlinear eqns
    [chis,~] = CorrelationEval(stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n),diffs,chis,params);
    ypfinal = phistobulkyseval(chis,stateVar(2*(n+1)+2+1:2*(n+1)+2+1+n),n,Vs,Pu,Pd,R,T,psat);
    resDAE(end-n:end-1) = stateVar(end-n:end-1)-ypfinal; %LFM eqns
    resDAE(end) = 1-sum(ypfinal); %LFM eqns

%Constant Perm (Usual Model)
    % res(1:2) = D(1:2)*W(5)+D(5)*W(1:2)+width*0.1.*(W(1:2)-W(6:7)); %feed species balance
    % res(3) = D(3)*W(5)+D(5)*W(3)+width*2*(W(3)-W(8)); %feed species balance
    % res(4) = D(4)*W(5)+D(5)*W(4);
    % res(5) = D(5)*sum(W(1:4))+width*sum(0.1*(W(1)-W(6))+0.1*(W(2)-W(7))+2*(W(3)-W(8))); %total feed balance
    % res(6:7) = D(6:7)*W(10)+D(10)*W(6:7)-width*0.1.*(W(1:2)-W(6:7)); %perm balance
    % res(8) = D(8)*W(10)+D(10)*W(8)-width*2*(W(3)-W(8)); %perm balance
    % res(9) = D(9)*W(10)+D(10)*W(9);
    % res(10) = D(10)*sum(W(6:9))-width*sum(0.1*(W(1)-W(6))+0.1*(W(2)-W(7))+2*(W(3)-W(8))); %total perm balance
end

