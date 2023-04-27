%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalMSA_RHS(phiFeed,localCompFlux,fugFeed,params)   %
% Description: RHS function for multiple shooting algorithm of asyMem     %
%                local flux.                                              %
% Input:       phiFeed       - n+1 dimensional vector of feed side        %
%                                membrane phase volume fractions          %
%              localCompFlux - n+1 dimensional vector of support layer    %
%                                compositions and total local mem flux    %
%              fugFeed       - n dimensional vector of penetrant          %
%                                feed side fugacities (torr)              %
%              params        - struct of system parameters                %
%                                (see dataBank function for specs)        %
% Output:      funRHS        - function value vector for nonlinear solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   [funRHS] = asyMemLocalMSA_RHS_pervap(phiFeed,localCompFlux,fugFeed,params,N)

%------------------------------------------------------------------------------------------------------------------------------------% 
%unpack parameters and modify as needed
    T = params.T;
    R = params.R;
    psat = params.psat;
    Vs = params.Vs;
    Pu = params.Pu;
    Pd = params.Pd;
    n = params.n;
    diffs = params.diffs;
    if params.memPhaseModel == 1 || params.memPhaseModel == 3
        chis = params.chis; %note chi_ji = chi_ij
    elseif params.memPhaseModel == 2
        chis = zeros(params.n+1);
    end

    if  params.swlDiffModel == 1
        params.Bffvtype = 'NA';
    elseif  any(params.Bffv == 0)
        params.Bffvtype = 'some0some1';
    else
        params.Bffvtype = 'cons1';
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%IVP solve
    shootingPoints = [];
    for i = 1:N
        zspan=[(i-1)*params.lmem/N params.lmem*(i/N)]; %um
        if params.memPhaseModel == 1
            diffVar = 0;
          % ode15s route
            funIVP = @(z,stateVar)asyMemLocal_IVP_pervap(z,stateVar,diffVar,localCompFlux,params); %solve IVP
            if i == 1
                [z,stateVars] = ode15s(funIVP,zspan,phiFeed);   
            else
                [z,stateVars] = ode15s(funIVP,zspan,...
                    localCompFlux(n+1+(i-2)*(n+1):(i-1)*(n+1)+n)); 
            end
          % manual discritization IVP solve
            %[z,stateVar] = asyMemLocalSA_IVP_TimeStep(phiFeed,localCompFlux,params);
            if i == N
                stateFinal = stateVars(end,1:n+1).';
            end
        elseif params.memPhaseModel == 2 || params.memPhaseModel == 3
          % DAE ode15i route  
            funIVP = @(z,W,D)asyMemLocal_IVP_pervap(z,W,D,localCompFlux,params); %solve IVP
            D0 = ones(n+1+n,1); %FH-DSM
            opt = odeset('InitialSlope', D0,'RelTol',1e-6,'AbsTol',1e-8);
             if i == 1
                [W0_new,D0_new] = decic(funIVP,0,[phiFeed;fugFeed],...
                    [ones(n+1,1);zeros(n,1)],D0,zeros(1,n+1+n),opt); %FH-DSM  
            else
                [W0_new,D0_new] = decic(funIVP,0,localCompFlux(n+1+(i-2)*(2*n+1):(i-1)*(2*n+1)+n),...
                    [ones(n+1,1);zeros(n,1)],D0,zeros(1,n+1+n),opt); %FH-DSM
             end
            opt = odeset(opt,'InitialSlope', D0_new,'RelTol',1e-6,'AbsTol',1e-8);
            [z,stateVars] = ode15i(funIVP,zspan,W0_new,D0_new,opt);
            if i == N
                stateFinal = stateVars(end,1:n).';
            end
        end
        if i ~= N
            shootingPoints = [shootingPoints;[stateVars(end,:).']];
        end
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%fsolve residual functions
    funRHS(1:n) = stateFinal(1:n);
    if N > 1 && params.memPhaseModel == 1
        funRHS(n+1:(N-1)*(n+1)+n) = localCompFlux(n+1:(N-1)*(n+1)+n).' - shootingPoints.';
    elseif N > 1 && (params.memPhaseModel == 2 || params.memPhaseModel == 3)
        funRHS(n+1:(N-1)*(2*n+1)+n) = localCompFlux(n+1:(N-1)*(2*n+1)+n) - shootingPoints;
    end
    
  %Manual discritization IVP solve mult shooting points
   % F(n+2:n+2+(params.N-1)*(n)-1) = comp_perm_flux(n+2:n+2+(params.N-1)*(n)-1)-W(floor(params.N/2),1:end-1).'; %only for single 
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%run table .dat file writing
% if norm(funRHS)^2 < 1E-5
%     norm(funRHS)^2;
% end
%     if norm(funRHS)^2<2E-9
%          partialFlux = localCompFlux(n+1)*localCompFlux(1:n).*params.Vs(1:n)/sum(localCompFlux(1:n).*params.Vs(1:n));
%          noww = num2str(now);
%          m_width = '';
%          for i=1:n
%              m_width = strcat(m_width,'%d ;');
%          end
%          fid = fopen(strcat('\Users\dweber31\Dropbox (GaTech)\Membrane_Stuff\Lively_Group_Contrib\MemSim-FFV-ErrorProp-newPIM\'...
%              ,params.Bffvtype,'BFFV_type',num2str(params.swlDiffModel),'diffswellmod',params.mixID,'_CrossDiffModel',num2str(params.diffModel),'_MemSorpModel'...
%              ,num2str(params.memPhaseModel),'FicksOG',num2str(params.memPhaseModel_FicksOG),'noThermoCoupling',num2str(params.noThermoCoupling),'_',noww,'.dat'),'w');
%          fprintf(fid,'Molar Composition in Feed\n');
%          fprintf(fid,strcat('%d\n'),params.yf);
%          fprintf(fid,'Partial Flux (LMH)\n');
%          fprintf(fid,strcat('%d\n'),partialFlux);
%          fprintf(fid,'Molar Composition in Permeate (Support Layer)\n');
%          fprintf(fid,strcat('%d\n'),localCompFlux(1:n));
%          fprintf(fid,'Total Flux (LMH)\n');
%          fprintf(fid,strcat('%d\n'),localCompFlux(n+1));
%          fprintf(fid,'Form of Membrane Phase Sorption model\n');
%          fprintf(fid,strcat('%s\n'),params.memPhaseModel);
%          fprintf(fid,'Diffusivities (cm2/s)\n');
%          fprintf(fid,strcat(m_width,'%d\n'),params.diffs.');
%          if params.memPhaseModel == 1 || params.memPhaseModel == 3
%              fprintf(fid,'Chis\n');
%              fprintf(fid,strcat(m_width,'%d\n'),params.chis.');
%          end
%          if params.memPhaseModel == 2 || params.memPhaseModel == 3
%              fprintf(fid,'b_s (atm^-1)\n');
%              fprintf(fid,strcat('%d\n'),params.bs);
%              fprintf(fid,'C_h (cm3 solv/cm3 sys)\n');
%              fprintf(fid,strcat('%d\n'),params.Ch);
%          end
%          if params.memPhaseModel == 2
%              fprintf(fid,'k_s (cc/cc torr^-1)\n');
%              fprintf(fid,strcat('%d\n'),params.ks);
%          end
%          fprintf(fid,'p_sat (atm)\n');
%          fprintf(fid,strcat('%d\n'),params.psat);
%          fprintf(fid,'Molar Volumes (cm3/mol)\n');
%          fprintf(fid,strcat('%d\n'),params.Vs);
%          fprintf(fid,'Temperature (K)\n');
%          fprintf(fid,strcat('%d\n'),params.T);
%          fprintf(fid,'Upstream Pressure (atm)\n');
%          fprintf(fid,strcat('%d\n'),params.Pu);
%          fprintf(fid,'Downstream Pressure (atm)\n');
%          fprintf(fid,strcat('%d\n'),params.Pd);
%          fprintf(fid,'Membrane Thickness (cm)\n');
%          fprintf(fid,strcat('%d\n'),params.lmem);
%          fprintf(fid,'Thickness Profile\n');
%          fprintf(fid,strcat(m_width,'%d\n'),z);
%          fprintf(fid,'Volume Fraction Profile\n');
%          fprintf(fid,strcat(m_width,'%d\n'),stateVars(:,1:n+1).');
%          if params.memPhaseModel == 2 || params.memPhaseModel == 3
%             fprintf(fid,'Fugacity Profile (torr)\n');
%             fprintf(fid,strcat(m_width,'%d\n'),stateVars(:,n+2:end).');
%          end
%          fprintf(fid,'Feed Fugacity (atm)\n');
%          fprintf(fid,'%d\n',fugFeed);
%          fprintf(fid,'B FFV\n');
%          fprintf(fid,strcat(m_width,'%d\n'),params.Bffv.');
%          if params.memPhaseModel == 2 || params.memPhaseModel == 3
%              fprintf(fid,'Permeate Fugacity (atm)\n');
%              fprintf(fid,'%d\n',fugFinal);
%          end
%          fclose(fid);
%     end
%------------------------------------------------------------------------------------------------------------------------------------% 

%% plots, activity, ect.
%      act_mem = [];
%      for i = 1:size(W,1)
%          phi_i = W(i,1:n+1).';
%          act_mem_i = [phis_to_bulk_y_s_eval(Chis,phi_i,n,V_s,P_u,P_d,R,T,p_sat).*...
%              exp(-V_s(1:n).*(P_u-2*P_d+p_sat)./(R*T))].';
%          act_mem = [act_mem;act_mem_i];
%      end
%final plot criteria
%       if norm(F)<1E-9
%           W_BVP = [0.0407532870000000,0.0204332710000000,0.00567018000000000,0.0278421050000000,0.0447313430000000,0.00440592700000000,0.000563365000000000,0.00342240000000000,0.000584862000000000,0.851593261000000;0.0407691770000000,0.0202406830000000,0.00551602300000000,0.0277219520000000,0.0445857100000000,0.00428170600000000,0.000548767000000000,0.00329775900000000,0.000564748000000000,0.852473470000000;0.0407866670000000,0.0200485250000000,0.00535991900000000,0.0276013730000000,0.0444410210000000,0.00415643500000000,0.000534098000000000,0.00317051700000000,0.000544174000000000,0.853357276000000;0.0408041520000000,0.0198563630000000,0.00520381700000000,0.0274807940000000,0.0442963240000000,0.00403116500000000,0.000519428000000000,0.00304328100000000,0.000523597000000000,0.854241071000000;0.0408232560000000,0.0196646440000000,0.00504574800000000,0.0273597890000000,0.0441525940000000,0.00390483500000000,0.000504686000000000,0.00291340000000000,0.000502553000000000,0.855128502000000;0.0408423550000000,0.0194729200000000,0.00488768000000000,0.0272387830000000,0.0440088570000000,0.00377850300000000,0.000489943000000000,0.00278352200000000,0.000481503000000000,0.856015921000000;0.0408630920000000,0.0192816540000000,0.00472762400000000,0.0271173520000000,0.0438661100000000,0.00365110200000000,0.000475129000000000,0.00265095300000000,0.000459978000000000,0.856907016000000;0.0408838250000000,0.0190903830000000,0.00456756800000000,0.0269959180000000,0.0437233550000000,0.00352369900000000,0.000460312000000000,0.00251838300000000,0.000438446000000000,0.857798099000000;0.0409062130000000,0.0188995820000000,0.00440550200000000,0.0268740590000000,0.0435816140000000,0.00339521500000000,0.000445424000000000,0.00238307600000000,0.000416430000000000,0.858692897000000;0.0409285980000000,0.0187087780000000,0.00424343500000000,0.0267521960000000,0.0434398650000000,0.00326672900000000,0.000430534000000000,0.00224776400000000,0.000394404000000000,0.859587683000000;0.0409526570000000,0.0185184580000000,0.00407933700000000,0.0266299090000000,0.0432991530000000,0.00313715200000000,0.000415572000000000,0.00210966600000000,0.000371886000000000,0.860486224000000;0.0409767140000000,0.0183281330000000,0.00391523600000000,0.0265076170000000,0.0431584340000000,0.00300757000000000,0.000400607000000000,0.00197156100000000,0.000349357000000000,0.861384755000000;0.0410024650000000,0.0181383080000000,0.00374908200000000,0.0263848990000000,0.0430187760000000,0.00287688800000000,0.000385571000000000,0.00183062000000000,0.000326325000000000,0.862287079000000;0.0410282130000000,0.0179484780000000,0.00358292600000000,0.0262621780000000,0.0428791110000000,0.00274620100000000,0.000370532000000000,0.00168967000000000,0.000303282000000000,0.863189394000000;0.0410556730000000,0.0177591610000000,0.00341469300000000,0.0261390290000000,0.0427405310000000,0.00261440200000000,0.000355421000000000,0.00154583300000000,0.000279725000000000,0.864095544000000;0.0410831330000000,0.0175698410000000,0.00324645700000000,0.0260158770000000,0.0426019460000000,0.00248259800000000,0.000340307000000000,0.00140198500000000,0.000256156000000000,0.865001685000000;0.0411123240000000,0.0173810470000000,0.00307612200000000,0.0258922970000000,0.0424644690000000,0.00234967000000000,0.000325120000000000,0.00125519700000000,0.000232063000000000,0.865911702000000;0.0411415150000000,0.0171922520000000,0.00290578300000000,0.0257687130000000,0.0423269890000000,0.00221673700000000,0.000309931000000000,0.00110839700000000,0.000207958000000000,0.866821712000000;0.0411724570000000,0.0170039960000000,0.00273332300000000,0.0256447010000000,0.0421906400000000,0.00208266800000000,0.000294668000000000,0.000958604000000000,0.000183316000000000,0.867735638000000;0.0412033990000000,0.0168157400000000,0.00256085700000000,0.0255206850000000,0.0420542900000000,0.00194859400000000,0.000279403000000000,0.000808796000000000,0.000158663000000000,0.868649561000000;0.0412361120000000,0.0166280380000000,0.00238624600000000,0.0253962400000000,0.0419190960000000,0.00181337200000000,0.000264064000000000,0.000655941000000000,0.000133460000000000,0.869567439000000;0.0412688270000000,0.0164403370000000,0.00221162900000000,0.0252717910000000,0.0417839020000000,0.00167814600000000,0.000248723000000000,0.000503070000000000,0.000108248000000000,0.870485318000000;0.0413033310000000,0.0162532040000000,0.00203484400000000,0.0251469120000000,0.0416498880000000,0.00154175900000000,0.000233307000000000,0.000347094000000000,8.25000000000000e-05,0.871407192000000;0.0413378390000000,0.0160660730000000,0.00185805300000000,0.0250220300000000,0.0415158780000000,0.00140536800000000,0.000217890000000000,0.000191104000000000,5.67000000000000e-05,0.872329072000000;0.0413741570000000,0.0158795240000000,0.00167906800000000,0.0248967160000000,0.0413830720000000,0.00126780400000000,0.000202396000000000,3.19000000000000e-05,3.03000000000000e-05,0.873254988000000];
%           W_SA = [0.0407532870000000,0.0204332710000000,0.00567018000000000,0.0278421050000000,0.0447313430000000,0.00440592700000000,0.000563365000000000,0.00342240000000000,0.000584862000000000,0.851593261000000;0.0407849480000000,0.0200781870000000,0.00538263800000000,0.0276195700000000,0.0444637290000000,0.00417497400000000,0.000536299000000000,0.00318841700000000,0.000547041000000000,0.853224198000000;0.0408192970000000,0.0197238360000000,0.00509182200000000,0.0273963250000000,0.0441977140000000,0.00394225500000000,0.000509113000000000,0.00295003400000000,0.000508435000000000,0.854861171000000;0.0408564190000000,0.0193702520000000,0.00479763000000000,0.0271723540000000,0.0439333670000000,0.00370771600000000,0.000481802000000000,0.00270708400000000,0.000469010000000000,0.856504366000000;0.0409081640000000,0.0189117890000000,0.00441108600000000,0.0268804130000000,0.0435920020000000,0.00340065600000000,0.000446158000000000,0.00238558000000000,0.000416739000000000,0.858647413000000;0.0409643560000000,0.0184545970000000,0.00401913300000000,0.0265873280000000,0.0432533690000000,0.00309068400000000,0.000410316000000000,0.00205666300000000,0.000363132000000000,0.860800422000000;0.0410252980000000,0.0179987890000000,0.00362140400000000,0.0262930320000000,0.0429176890000000,0.00277760400000000,0.000374265000000000,0.00171976800000000,0.000308085000000000,0.862964067000000;0.0410911450000000,0.0175444380000000,0.00321771400000000,0.0259975000000000,0.0425850980000000,0.00246132300000000,0.000337998000000000,0.00137457100000000,0.000251534000000000,0.865138679000000;0.0411619830000000,0.0170916000000000,0.00280796100000000,0.0257007260000000,0.0422556940000000,0.00214178600000000,0.000301515000000000,0.00102085700000000,0.000193433000000000,0.867324446000000;0.0412378840000000,0.0166403280000000,0.00239206200000000,0.0254027060000000,0.0419295670000000,0.00181895100000000,0.000264812000000000,0.000658430000000000,0.000133741000000000,0.869521518000000;0.0413189120000000,0.0161906700000000,0.00196994000000000,0.0251034380000000,0.0416067980000000,0.00149278300000000,0.000227890000000000,0.000287115000000000,7.24000000000000e-05,0.871730032000000];
%      
%          z_SA = linspace(0,params.lmem,size(W_SA,1));
%          z_BVP = linspace(0,params.lmem,size(W_BVP,1));
%          P = W(:,1);
%          %P_lin = [phi_s(1,1),phi_s(end,1)];
%          P1_SA = W_SA(:,1);
% %         %P1_lin = [phi_s(1,2),phi_s(end,2)];
%          P2_SA = W_SA(:,2);
% %         %P2_lin = [phi_s(1,3),phi_s(end,3)]';
%          P3_SA = W_SA(:,3);
%          P4_SA = W_SA(:,4);
%          P5_SA = W_SA(:,5);
%          P6_SA = W_SA(:,6);
%          P7_SA = W_SA(:,7);
%          P8_SA = W_SA(:,8);
%          P9_SA = W_SA(:,9);
%          P10_SA = W_SA(:,10)
%          P1_BVP = W_BVP(:,1);
% %         %P1_lin = [phi_s(1,2),phi_s(end,2)];
%          P2_BVP = W_BVP(:,2);
% %         %P2_lin = [phi_s(1,3),phi_s(end,3)]';
%          P3_BVP = W_BVP(:,3);
%          P4_BVP = W_BVP(:,4);
%          P5_BVP = W_BVP(:,5);
%          P6_BVP = W_BVP(:,6);
%          P7_BVP = W_BVP(:,7);
%          P8_BVP = W_BVP(:,8);
%          P9_BVP = W_BVP(:,9);
%          P10_BVP = W_BVP(:,10)
% %         
%          figure(2)
%          %plot(z,P,"-o",z_lin,P_lin);%,z,P1,"-s",z_lin,P1_lin,z,P2,"-*",z_lin,P2_lin
%          plot(z_SA,P1_SA,'o-',z_BVP,P1_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(3)
%          plot(z_SA,P2_SA,'o-',z_BVP,P2_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(4)
%          plot(z_SA,P3_SA,'o-',z_BVP,P3_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(5)
%          plot(z_SA,P4_SA,'o-',z_BVP,P4_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(6)
%          plot(z_SA,P5_SA,'o-',z_BVP,P5_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(7)
%          plot(z_SA,P6_SA,'o-',z_BVP,P6_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(8)
%          plot(z_SA,P7_SA,'o-',z_BVP,P7_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(9)
%          plot(z_SA,P8_SA,'o-',z_BVP,P8_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(10)
%          plot(z_SA,P9_SA,'o-',z_BVP,P9_BVP);
%                   title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
%          figure(11)
%          plot(z_SA,P10_SA,'o-',z_BVP,P10_BVP);
%          title("MS Simulated-Volume Fraction Profile");
%          xlabel("z (cm)");
%          ylabel("Volume Faction of Comp_i");
%          legend("SA","BVP");
% %         %legend("MS Simulation-Toluene","MS Simulation-iso-Octane","MS Simulation-iso-cetane");,"MS Simulation-n-Octane","Line","MS Simulation-1-MN","Line"
%      end



    



    
    
   
    
    
    
    
    

    
    
    
    
     

    
    
    

    
    

