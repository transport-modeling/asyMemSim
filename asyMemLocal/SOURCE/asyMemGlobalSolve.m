%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSolve(sysInfoExt or n/a)                        %
% Description: Solve global flat plat cocurrent asymetric membrane module %
%                model problem.                                           %
% Input:       params  - struct of system parameters                      %
%                          (see dataBank function for specs)              %
% Output:      flowVec - vector containg final component flow rates (L/h) %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [flowVec] = asyMemGlobalSolve(params)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters
    n = params.n;
    xSpan = [0 1]; %m
    cFeedInt = [params.y_f;0.01];
    vTFeedInt = 10; %L/h
    cPermInt = [zeros(n,1);1];
    vTPermInt = 10; %L/h
    phiPermGuess = [params.y_f*0.01;1-sum(params.y_f*0.01)];
    yPermGuess = cFeedInt(1:n)/sum(cFeedInt);
    nTotGuess = 10; %LMH
    dcFeedGuess = [-0.01*ones(n,1);0.01];
    dvTFeedGuess = -ones(1,1);
    dcPermGuess = [0.01*ones(n,1);-0.01];
    dvTPermGuess = ones(1,1);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%solve DAE
    fun = @(x,W,D)asyMemGlobalRHS(x,W,D,params); 
    % D0 = [dcFeedGuess;dvTFeedGuess;dcPermGuess;dvTPermGuess]; 
    % W0 = [cFeedInt;vTFeedInt;cPermInt;vTPermInt];

    D0 = [dcFeedGuess;dvTFeedGuess;dcPermGuess;dvTPermGuess;zeros(2*(n+1),1)]; 
    W0 = [cFeedInt;vTFeedInt;cPermInt;vTPermInt;phiPermGuess;...
        yPermGuess;nTotGuess];
    opt = odeset('InitialSlope', D0);
    [W0_new,D0_new] = decic(fun,0,W0,[ones(2*(n+1)+2,1);zeros(2*(n+1),1)],D0,zeros(4*(n+1)+2,1),opt);
    opt = odeset(opt,'InitialSlope', D0_new);
    [x,W]=ode15i(fun,xSpan,W0_new,D0_new,opt);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%9-comp plots
    figure(1)
    plot(x,W(:,1)./sum(W(:,1:n+1),2),"-o",...
        x,W(:,2)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,3)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,4)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,5)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,6)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,7)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,8)./sum(W(:,1:n+1),2),"-o"...
        ,x,W(:,9)./sum(W(:,1:n+1),2),"-o");
    title("MS Global Module Feed Channel Compositions");
    axis([0 1 0 0.3]);
    xlabel("z (m)");
    ylabel("Composition (mol/mol)");
    
%     figure(2)
%     plot(z,W(:,1).*W(:,11),"-o",z,W(:,2).*W(:,11),z,W(:,3).*W(:,11),z,W(:,4).*W(:,11),"-o",z,W(:,5)...
%         .*W(:,11),z,W(:,6).*W(:,11),z,W(:,7).*W(:,11),"-o",z,W(:,8).*W(:,11),z,W(:,9).*W(:,11));
%     title("MS Global Module Feed Channel Component Flow Rates");
%     xlabel("z (m)");
%     ylabel("Component Flow Rate (mol/hr)");
    feedTotFlow_fin = W(end,1).*W(end,11)+W(end,2).*W(end,11)+W(end,3).*W(end,11)+...
        W(end,4).*W(end,11)+W(end,5).*W(end,11)+W(end,6).*W(end,11)+...
        W(end,7).*W(end,11)+W(end,8).*W(end,11)+W(end,9).*W(end,11)+W(end,10).*W(end,11)
%     
%     figure(3)
%     plot(z,W(:,12).*W(:,22),"-o",z,W(:,13).*W(:,22),z,W(:,14).*W(:,22),z,W(:,15).*W(:,22),"-o",z,W(:,16)...
%         .*W(:,22),z,W(:,17).*W(:,22),z,W(:,18).*W(:,22),"-o",z,W(:,19).*W(:,22),z,W(:,20).*W(:,22));
%     title("MS Global Module Peremate Channel Component Flow Rates");
%     xlabel("z (m)");
%     ylabel("Component Flow Rate (mol/hr)");
%     permTotFlow = W(end,12).*W(end,22)+W(end,13).*W(end,22)+W(end,14).*W(end,22)+...
%         W(end,15).*W(end,22)+W(end,16).*W(end,22)+W(end,17).*W(end,22)+...
%         W(end,18).*W(end,22)+W(end,19).*W(end,22)+W(end,20).*W(end,22)+W(end,21).*W(end,22)
    
    figure(4)
    plot(x,W(:,end).*W(:,end-9).*params.V_s(1)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-8).*params.V_s(2)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-7).*params.V_s(3)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-6).*params.V_s(4)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-5).*params.V_s(5)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-4).*params.V_s(6)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-3).*params.V_s(7)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-2).*params.V_s(8)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o",x,...
        W(:,end).*W(:,end-1).*params.V_s(9)./sum(W(:,end-9:end-1).*(params.V_s(1:n)).',2),"-o");
    axis([0 1 0 0.55]);
    title("MS Global Module Partial Transmembrane Fluxes");
    xlabel("z (m)");
    ylabel("Component Partial Flux (LMH)");   
    totFluxFin = W(end,end)
    totFLuxInt = W(1,end)
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
    %3-comp plots
%     figure(1)
%     plot(z,W(:,1)./sum(W(:,1:n+1),2),"-x",z,W(:,2)./sum(W(:,1:n+1),2),"-x",z,W(:,3)./sum(W(:,1:n+1),2),"-x");
%     title("MS Global Module Feed Channel Compositions");
%         axis([0 1 0.24 0.4]);
%     xlabel("z (m)");
%     ylabel("Composition (mol/mol)");
% 
% %     figure(2)
% %     plot(z,W(:,1).*W(:,5),"-o",z,W(:,2).*W(:,5),"-o",z,W(:,3).*W(:,5),"-o");%,z,P1,"-s",z_lin,P1_lin,z,P2,"-*",z_lin,P2_lin
% %     title("MS Global Module Feed Channel Component Flow Rates");
% %     xlabel("z (m)");
% %     ylabel("Component Flow Rate (mol/hr)");
%      feedTotFlow_fin = W(end,1).*W(end,5)+W(end,2).*W(end,5)+W(end,3).*W(end,5)+W(end,4).*W(end,5)
% %     
% %     figure(3)
% %     plot(z,W(:,6).*W(:,10),"-o",z,W(:,7).*W(:,10),"-o",z,W(:,8).*W(:,10),"-o");
% %     title("MS Global Module Permeate Channel Component Flow Rates");
% %     xlabel("z (m)");
% %     ylabel("Component Flow Rate (mol/hr)");
% %     permTotFlow = W(end,6).*W(end,10)+W(end,7).*W(end,10)+W(end,8).*W(end,10)+W(end,9).*W(end,10)
% %     
%     figure(4)
%     plot(z,W(:,end).*W(:,end-3).*params.V_s(1)./sum(W(:,end-3:end-1).*(params.V_s(1:n)).',2),"-x",z,...
%         W(:,end).*W(:,end-2).*params.V_s(2)./sum(W(:,end-3:end-1).*(params.V_s(1:n)).',2),"-x",z,...
%         W(:,end).*W(:,end-1).*params.V_s(3)./sum(W(:,end-3:end-1).*(params.V_s(1:n)).',2),"-x");
%     axis([0 1 0.04 0.2]);
%     title("MS Global Module Partial Transmembrane Fluxes");
%     xlabel("z (m)");
%     ylabel("Component Partial Flux (LMH)");
%     totFluxFin = W(end,end)
%     totFLuxInt = W(1,end)
%------------------------------------------------------------------------------------------------------------------------------------% 

end

