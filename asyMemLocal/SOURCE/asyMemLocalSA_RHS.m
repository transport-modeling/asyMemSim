%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemLocalSA_RHS(phiFeed,localCompFlux,fugFeed,params)    %
% Description: RHS function for shooting algorithm of asyMem local flux.  %
% Input:       phiFeed       - n+1 dimensional vector of feed side        %
%                                      membrane phase volume fractions    %
%              localCompFlux - n+1 dimensional vector of support layer    %
%                                compositions and total local mem flux    %
%              fugFeed       - n dimensional vector of penetrant          %
%                                     feed side fugacities (torr)         %
%              params        - struct of system parameters                %
%                                     (see dataBank function for specs)   %
% Output:      funRHS        - function value vector for nonlinear solver %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function   [funRHS] = asyMemLocalSA_RHS(phiFeed,localCompFlux,fugFeed,params)

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
    zspan=[0 params.lmem]; %um
    if params.memPhaseModel == 1 && params.noGammaFugacityODEs == 0
      % ode15s route
        diffVar = 0;
        funIVP = @(z,stateVar)asyMemLocal_IVP(z,stateVar,diffVar,localCompFlux,params); %solve IVP
        opt = odeset('RelTol',1e-7,'AbsTol',1e-9);
        [z,stateVars] = ode15s(funIVP,zspan,phiFeed,opt);
%         [z,stateVars] = ode15s(funIVP,zspan,1000*phiFeed,opt);
      % manual discritization IVP solve
        %[z,stateVar] = asyMemLocalSA_IVP_TimeStep(phiFeed,localCompFlux,params);
      % DAE route (overkill)...also need to change IVP function
%         funIVP = @(z,W,D)asyMemLocal_IVP(z,W,D,localCompFlux,params); %solve IVP
%         D0 = ones(n+1,1); %FH-DSM
%         opt = odeset('InitialSlope', D0,'RelTol',1e-4,'AbsTol',1e-6);
%         [W0_new,D0_new] = decic(funIVP,0,[phiFeed],[ones(n+1,1)],D0,zeros(1,n+1),opt); %FH-DSM
%         opt = odeset(opt,'InitialSlope', D0_new,'RelTol',1e-4,'AbsTol',1e-6);
%         [z,stateVars] = ode15i(funIVP,zspan,W0_new,D0_new,opt);   
        phiFinal = stateVars(end,1:n+1).';
%         phiFinal = stateVars(end,1:n+1).'/1000;
        [params.chis,~] = correlationEval(phiFinal,diffs,chis,params);
        ypFinal = phis2yPhaseEq_FH_RHS(phiFinal,params);% FH based
    elseif params.memPhaseModel == 2 || params.memPhaseModel == 3 || (params.memPhaseModel == 1 && params.noGammaFugacityODEs == 1)
      % DAE ode15i route  
        funIVP = @(z,W,D)asyMemLocal_IVP(z,W,D,localCompFlux,params); %solve IVP
        D0 = ones(n+1+n,1); %FH-DSM
%         D0 = [-371;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;-6.81267603133497e-07;-1.66766081827139e-05;-1.97811174016514e-06;-1.06423157442261e-05;-4.49217506089836e-06;-0.000236446622315984;-3.11956467065235e-05;-3.04193368480250e-05;-5.92319234121443e-06;-0.000115124094765646;-7.70511746214574e-08;-6.56447161078459e-05;4.94643875974396e-10;-0.000178422301792449;-0.00162704162933389;1.75134369582768e-08;-0.000208516460769752;-5.43465458889041e-08;-1.18254693395497e-05;1.17425888505853e-08;-4.17980277022177e-06;-5.02090160165677e-08;-6.46216940232508e-07;7.66567980813424e-09;-3.12568445440392e-06;-4.03618617312210e-06;-8.89118902506012e-06;-0.000557517870519425;-0.00110447160762660;-4.82675387094053e-07;-3.57246617847132e-06;-8.23342422699900e-05;-0.000143646491973842;2.67449402656402e-08;-2.21904439032575e-05;-4.36165878104813e-05;-1.34487988504194e-05;-0.000239397873041127;-3.40039427255860e-07;-2.91945846164660e-05;-1.80336910723192e-05;-1.93347347639739e-08;-0.000244506844241770;-4.62703834994047e-05;-2.25026195000286e-05;-1.56074097514411e-05;-2.00479156234343e-07;1.57621190762970e-08;-0.000744790973056243;-0.000174680598497372;-0.000883401244865736;-7.34910149608115e-06;-3.51956037687479e-05;-5.81512082710015e-06;-7.39340765212410e-07;-9.58491300909072e-05;-3.25474535253315e-05;-0.000454508175117547;-2.59641027331662e-05;-4.93040829452616e-05;1.17340941908188e-09;-0.000663124140098064;-7.01137845493142e-05;-1.01845246912062e-05;-1.84881057606083e-06;-0.000397441227646769;-0.000352889220730773;-9.41340296286461e-05;-0.000248473732409956;-8.69526002879305e-07;-2.10430399328313e-05;-4.72984787360007e-07;-0.000208223693379007;-6.83281881420481e-06;-0.000181761127032851;-3.07659961865815e-08;-9.41509126251062e-05;-5.53112290447502e-05;-7.69373480253084e-06;-4.14311282554797e-06;-3.02645625225828e-05;-5.29329531725070e-08;-0.000613013148394498;-5.29891502657863e-05;-3.05799405274549e-06;-0.000112270742503507;-0.000546259513798992;-4.05284765549250e-07;1.42615693204600e-08;-3.58087077395157e-07;-1.46901953128597e-05;-1.02516311568746e-06;-2.10832875461854e-05;-2.07953534055980e-07;-2.66856522816007e-05;-2.13929892326208e-06;-0.000143783176271631;-0.000263210299571189;-8.25197606719079e-07;-4.72845684219964e-06;-3.87079687209214e-05;-2.64095312362457e-06;-1.25355233276932e-06;-6.08707462788093e-05;-6.71590053774899e-05;-1.05445316414845e-08;-2.39563337567447e-07;-5.70536731743799e-06;-2.44910211749946e-06;-2.52484559566329e-05;-2.32221393040894e-05;-0.000156442375964759;-0.000156773047348708;-8.67704783376042e-06;-0.000436450496183670;-8.86764945331217e-06;-1.47592594141240e-06;3.11991265482429e-09;1.37657281086948e-08;-0.000442226243992553;1.63974891137447e-08;-4.59806546195201e-06;-6.73105654069167e-05;-1.29450646788935e-05;-0.000957308263655893;-1.57856714519415e-05;-8.08967046155059e-06;-9.92931478643660e-07;-5.55889638958949e-07;-6.47521424984886e-06;-0.000595162450981123;-3.87627508882445e-05;3.67887453650742e-09;-0.00139623640797349;-9.35961010539155e-07;-1.11080498972555e-06;-0.000492388218346711;-1.37448810441670e-07;-0.000201862357751286;-2.95050513849524e-06;-1.40380996126857e-08;-0.000365714030193345;-0.000165649740728034;9.58336631307928e-09;-6.05398920616093e-05;-7.73955796719403e-05;2.64860427071166e-09;-1.56367085981820e-08;-2.11129256820614e-05;-0.000206202610418962;-2.98477357876068e-05;-0.000659812046941656;-0.000322896942389543;-7.17512095042165e-08;-5.66315141490168e-06;-1.69128125525407e-07;-7.74311218260108e-05;-0.000227609346756196;-0.000148852697109753;-9.94658044327925e-06;-4.07068879183438e-05;-2.23840759585375e-05;-1.94572283179851e-05;-5.35143740850222e-05;-6.79428144540661e-05;-0.000349298456218072;-6.43811899456654e-06;2.95832419439232e-08;-0.000107237825882150;-8.45659299106297e-06;-8.07344229359296e-05;-2.08633986563568e-06;-6.54100748563149e-06;-2.38818585372731e-08;-4.00545291609764e-08;2.41281382990522e-08;-0.000135499977580619;1.10509477308823e-08;-0.000613973278135584;-7.95396862335574e-05;-5.52944594387895e-06;-4.03352029612164e-05;-8.86412815901290e-06;-1.85633406495913e-06;-1.89549357824720e-05;9.15717166725987e-09;-0.000104385448235515;1.05011195355615e-11;1.65639322297433e-08;-7.58315701990119e-06;-2.39172481228793e-06;-1.08817256018486e-06;1.48452768688937e-08;-1.56974698904373e-05;-1.01536824619442e-07;1.62176662540516e-08;-5.13446315392905e-06;-1.10379572752810e-05;-0.000456389128211963;2.54461137217145e-08;-0.000486474763406139;2.24512497209316e-08;-0.000500853203566844;-0.000434369507762220;-2.84338570744545e-05;-8.05171174392382e-05;1.16185705156451e-08;-5.17637230319016e-07;3.61125312647838e-09;-5.90882239665468e-07;-9.27371437700436e-05;-0.000124606756999032;-2.00081364564635e-07;-0.000250686720890812;5.34533890753743e-09;-3.33475317401417e-05;-5.84228570244017e-05;-1.55970682393657e-07;-1.84251924854141e-07;-0.000194606684654135;-7.82505977021288e-06;-0.000246615928826829;2.06950582705369e-08;-0.000223542934997000;-6.12788774649146e-05;-8.76761127959069e-05;-1.55377304748585e-05;-7.56356206737174e-05;-3.81462888092277e-05;-8.73938164927014e-05;-0.000174559624601875;-0.000160830684809816;-0.000113150352464198;-5.32622871920434e-05;-1.44280036041514e-06;-7.18240148952566e-08;-6.61886255690113e-05;-3.55597596822132e-07;-5.76085985270105e-05;-4.18161190635606e-06;-2.38784676527296e-07;-1.75577092061788e-05;-7.76541682019603e-06;-1.09162513501482e-05;-0.000349058562962849;-3.84845479549174e-05;-1.40288902093562e-05;-0.00168781796584200;-7.19474625325265e-06;2.22995155272688e-08;-4.31315690024652e-05;-0.000134050165770458;-6.99586555545188e-06;-0.000112469076267803;-0.00249898546584206;-0.000303971874629035;1.89866836766451e-08;-0.000924148055885934;-0.00226588417012337;-0.000123683612670875;-1.35438184441353e-06;-1.65409126818538e-06;-3.71842851178982e-06;-1.81854378403778e-05;-3.59763635548336e-05;-0.000437481379945506;-0.000382885424800237;-2.36149389448661e-05;-2.86112013072506e-08;-4.77920316035129e-05;-1.33075106489222e-05;-0.000154550774470951;-3.41710444923991e-06;-2.82550026858943e-07;-3.79101965681902e-05;-3.04069836591744e-05;1.42222372852778e-08;-6.97647299423203e-09;-8.89704872237073e-05;-0.000182249398682191;1.00812426309150e-08;-2.98855571297583e-06;-0.000329228675691723;1.77423427925493e-08;-4.30467249110885e-06;-2.17442913115568e-06;-1.67470572544800e-05;-7.58929873548989e-06;1.15942132149363e-08;-5.67115212053495e-06;-0.000474971769989839;-1.95619032062556e-06;-2.20122702496898e-06;-2.42623823064828e-06;2.27280283038307e-08;2.54669591959563e-08;-0.000138583265346788;-0.000112123609372073;-3.53440704600908e-05;-0.000196766569416456;-6.55910487204723e-05;-2.06995071070801e-05;-5.07070199390102e-06;-0.00163751364429364;-7.30481293971841e-05;-3.00846653033992e-05;-0.000293639928311633;-2.70527515796063e-07;-2.60266842342334e-05;-4.26159869217189e-06;-0.000654790871023842;2.01499152979624e-08;-1.05407169366070e-05;-1.80769779218590e-06;-2.20261485715365e-06;-3.10090890832370e-06;-9.28634422539518e-06;-7.54903041153844e-08;-1.03747715540658e-10;-1.49653917699522e-05;1.53837984982924e-08;-5.99999314541299e-06;1.60996958735097e-08;-9.24853900716081e-06;-2.68985100580566e-05;-0.000588241136224488;-4.17709575547540e-06;-1.26800746536388e-05;-2.68487991711079e-06;-3.59481184024741e-05;-5.67721003059128e-09;-0.000205156723665223;-1.08167313781047e-07;-1.03974498649139e-05;-0.00453020356744181;-6.14373564188363e-06;-0.000159108021547022;-0.000137358993933602;-9.82478073840382e-07;-0.000309858284877751;-1.42734920549548e-07;-0.00390967397167783;-3.13405118340281e-05;-5.04617681329692e-05;-1.20497828997660e-06;-2.99403763198111e-07;-4.77930086433498e-07;-6.09025018004594e-06;-9.42253685632261e-05;-0.000300719893604792;-1.09972460867769e-05;-3.49918115780847e-05;-1.20208503993025e-05;-3.91823913211766e-06;1.79989028618039e-08;-9.39758357014642e-06;-7.67404053268987e-07;-1.36183955713769e-05;-2.95318419592462e-05;-4.86432202422568e-05;-0.000133445248083170;-1.27927062041095e-07;-2.69807214714392e-08;-0.000585515627975635;-0.000128037438313269;-4.51682943939757e-05;-0.00477283661768443;-1.19397578928104e-06;-0.000157084671840480;-0.000427842420452931;1.69934892742750e-08];
        opt = odeset('InitialSlope', D0,'RelTol',1e-6,'AbsTol',1e-8);
        [W0_new,D0_new] = decic(funIVP,0,[phiFeed;fugFeed],[ones(n+1,1);zeros(n,1)],D0,zeros(1,n+1+n),opt); %FH-DSM
        opt = odeset(opt,'InitialSlope', D0_new,'RelTol',1e-6,'AbsTol',1e-8);
        [z,stateVars] = ode15i(funIVP,zspan,W0_new,D0_new,opt);
        fugFinal = (stateVars(end,n+2:end).');
%         fugFinal = (exp(stateVars(end,n+2:end)).');
         ypFinal = fugFinal./(exp(-Vs(1:n).*(Pu-Pd)/(R*T)))./params.purFug;
%        ypFinal = fugFinal./(exp(-Vs(1:n).*(Pu-Pd)/(R*T))); %if assuming fs = activity
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%fsolve residual functions
    funRHS(1:n) = localCompFlux(1:n)-ypFinal(1:n);
    funRHS(n+1) = 1-sum(ypFinal(1:n)); 
    
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
%              ,params.Bffvtype,'BFFV_type',num2str(params.swlDiffModel),'diffswellmod',strjoin(params.mixID),'_CrossDiffModel',num2str(params.diffModel),'_MemSorpModel'...
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



    



    
    
   
    
    
    
    
    

    
    
    
    
     

    
    
    

    
    
