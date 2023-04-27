%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    dataBank(sysInfo)                                          %
% Description: Build parameter matricies based on system specifications.  %
% Input:       sysInfoExt - external struct defining simulation specs     %
%                             (memID, mixID, yf, n, lmem, Pu, Pd, T, R    %
%                              memPhaseModel, diffModel, swlDiffModel)    %
% Output:      params     - struct of system parameters                   %
%                             (compID, n, Vs, HanSolParam, psat, chis [FH %
%                             or FH-LM], diffs, Ch & bs [FH-LM or DSM],   %
%                             ks [DSM], unitActPhis, Bffv, and all        %
%                             fields listed in sysInfo struct)            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [params] = dataBank_EXAMPLE(sysInfo)

%------------------------------------------------------------------------------------------------------------------------------------% 
%full parameter/property sets
    if contains(sysInfo.memID,'SBAD1')
      % Lively/Ronita's 9-comp data TOL/MCH/1MN/DCN/NOC/IOC/TBB/TPB/ICE/SBAD-1
        params.lmem = 0.3; % thickness of active membrane layer um
        params.compID = struct('TOL', 1, 'MCH', 2, 'MNP', 3, 'DEC', 4, 'NOC', 5, 'IOC', 6, 'TBB', 7, 'TPB', 8, 'ICT', 9);
        params.n = 9; 
        params.Vs = [106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;62326]; %cm3/mol
        params.HanSolParam = [18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0]; %[delD,delP,delH;...]
        params.psat = [28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458];%torr
        
        if contains(sysInfo.memPhaseModel,'F-H')    
            params.chis_im = [0.839;1.672;0.705;2.785;1.164;3.046;1.648;2.5;3.128];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [15.4638846409154;5.59021311138405;0.0535253755307862;1.66732636951057;5.74769726115217;6.9787;1.40290524341143;0.0692133718931809;0.589420725298558]; %um2/s FH Exp Based tweaked for single comp MS diffs 
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [14.1575747503318;2.87064997759363;0.0426220582870798;1.29680939844570;4.45897590219130;4.2057;0.975934082369876;0.0692133718931809;0.536239907969929]; %Mixed Exp (LOW ERROR) Based FH Tweaked to match single comp FH um2/s   
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0
                    params.diffs_im = [3.62215722661699;3.32944378753008;0.00962069984261326;1.39554404865361;2.36989760060502;6.09179492760562;0.833944291340556;0.0554924921379044;0.520245216257022]; %um2/s FICKS FH-BC tweaked for single comp 20bar
                elseif sysInfo.lowDiffErrorBar == 1    
                    params.diffs_im = [3.31617590838678;1.70971437749573;0.00766092765208338;1.08542314895282;1.83853042563447;3.67121416166937;0.580135159196378;0.0554924921379044;0.473305798253372]; %um2/s FICKS (LOW ERROR) FH-BC tweaked for single comp 20 bar
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.110;0.0996;0;0.0620;0;0.07;0.00975;0.0150;0.00998];
            params.bs = [0.810;0.11;212;0;0.06;0.68;0;0;46.5]; %torr^-1 mixed exp and genome
            params.ks = [0.0208;0.00079;19.3;0.00802;0.00987;0.00026;0.0568;1.13;0.251]; 
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [12.0068169139374;5.05131456788025;0.0343697334964408;5.16294418206992;9.36771416660342;1.5615;1.23897800067487;0.0548458701450260;0.485023478949192]; %Mixed Exp Based DSM Tweaked to match single comp backup1.22 i-OCT for 20 bar fit but 2.852 for fit of 30bar-60bar %um2/s
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [10.9925424245108;2.59391829377276;0.0273684914817629;4.01562325166401;7.26732982714958;0.9410;0.861897739597925;0.0548458701450260;0.441261962162586]; %Mixed Exp (LOW ERROR) Based DSM Tweaked to match single comp DSM um2/s
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0
                    params.diffs_im = [24.2137465214917;13.4274957298005;0.0760497427477021;5.20630692724413;10.8423340728343;8.71181973463473;1.39871601154965;0.0572833772672035;0.658451327584059]; %um2/s FICKS DSM-BC tweaked for single comp
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [22.1682930475097;6.89520050987365;0.0605581284842714;4.04934983120933;8.41131746831117;5.25016950894592;0.973019834118825;0.0572833772672035;0.599042185239266]; %um2/s FICKS (LOW ERROR) DSM-BC tweaked for single comp
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.chis_im = [1.01;2.32;0.79;3.06;1.25;3.76;1.69;2.66;3.63];
            params.Ch = [0.31;0.07302;0.481;0.00761;0.051;0.0112;0.00821;0.00689;0.00801];
            params.bs = [0.431;0.278;86.057;2.145;0.016;0.157;0.083;509.449;110.020]; %torr^-1 
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [11.4609361017534;5.30042942577299;0.0328169584688351;1.77959539744715;6.38183951490904;6.9251;1.47399503331522;0.0678674014003992;0.593101878582401]; %Mixed Exp Based FHLM Tweaked to match single comp FH-DSM um2/s
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [10.4927748326011;2.72184213749381;0.0261320224893682;1.38412975319562;4.95093379528772;4.1734;1.02538784926353;0.0678674014003992;0.539588927306691];
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [7.41480938962495;9.58113962154477;0.0149431737014035;1.84888715329978;2.98162315663263;12.6917932517404;0.902168297065590;0.0685073585122290;0.857099124490609]; %um2/s FICKS FHLM-BC tweaked for single comp
                elseif sysInfo.lowDiffErrorBar == 1 
                    params.diffs_im = [6.78844421295173;4.92004467148238;0.0118991938729840;1.43802334134041;2.31309778529407;7.64869659541521;0.627595337091547;0.0685073585122290;0.779766872654510]; %um2/s FICKS (LOW ERROR) FHLM-BC tweaked for single comp
                end
            end
        end
        
    %  FFV Diff Model
        params.unitActPhis = [0.420;0.1;0.55;0.0248;0.1569;0.02;0.0429;0.0288;0.0153]; %used with FFV calcs. Get these values from unit activity exp sorption volume fractions 
        params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1]; 

        
    elseif contains(sysInfo.memID,'PIM1')

      % Lively/Ronita's Exp Data for 5-comp TOL/HEP/PXY/OXY/ICE/PIM-1
        params.lmem = 1.5; %thickness of active membrane layer um
        params.compID = struct('TOL', 1, 'HEP', 2, 'PXY', 3, 'OXY', 4, 'ICT', 5);
        params.n = 5;
        params.Vs = [106.521;146.927;123.738;121.196;293.267;71429]; %cm3/mol
        params.HanSolParam = [18,1.4,2;15.3,0,0;17.6,1,3.1;17.8,1,3.1;16.3,0,0]; %[delD,delP,delH;...]
        params.psat = [28.998;44.854;8.803;6.637;0.0458]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [0.648;0.803;0.642;0.56;0.698];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [50.3125149850118;173.948724762266;35.0176589028337;20.2884438636580;2.79466382773699]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [37.3836039137959;141.294090512243;26.2347455779958;20.8208136270388;1.52857435515592]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [6.83432114205229;40.7634564837128;4.87536536909399;1.31261700013715;0.613893311027835]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [5.07809120717587;33.1111108621231;3.65033701457585;1.07169362670079;0.334010769242320]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.221;0.579;0.321;0;0.748];
            params.bs = [6880;0.973;15800;31.5;255]; %torr^-1 mixed exp and genome
            params.ks = [0.0525;0.00744;0.171;0.538;10.7];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.9500096480641;109.920896713247;17.9512856335608;5.59517317369966;1.36250783890009]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.7676366509148;89.2859269313939;13.4406834054801;4.56821100305585;0.741321469339973]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [89.6901094516630;531.240864816618;62.8682639525013;26.6727574346793;6.35724827979464]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [66.6422526266453;431.513339457563;47.0714159032490;21.7771246968768;3.45888992433533]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.Ch = [0.77;0.68;0.73;0.71;1.01];
            params.bs = [0.59;0.68;5.32;6.55;161.73];
            params.chis_im = [0.73;1.38;0.71;0.57;1.02];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [24.4320861716655;106.715078468667;16.5111963633180;6.11999794342432;1.16129980370385]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.1537212921656;86.6819229437563;12.3624442465211;6.11968813199987;0.631846701207788]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.5760737647037;278.508561986389;7.59319690360824;1.54591628595111;2.06900952589494]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.60134564606960;226.225366152481;5.68526164101519;1.26217215750965;1.12571916140146]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.65;0.46;0.67;0.85;0.6];
    params.Bffv = 0.03*[1;1;1;1;1]; 
    
    elseif contains(sysInfo.memID,'Albatross-5-12comp')

        params.lmem = 0.3; %thickness of active membrane layer um
        %params.compID = struct('TOL', 1, 'MCH', 2, 'MNP', 3, 'DEC', 4, 'NOC', 5, 'IOC', 6, 'TBB', 7, 'TPB', 8, 'ICT', 9);
        params.compID = struct('IOC', 1, 'TOL', 2, 'BCH', 3, 'TET', 4, 'MNP', 5, 'ICT', 6, 'TPB', 7, 'PRI', 8, 'DDB', 9, 'PYR', 10, 'DCS', 11, 'SQU', 12);
        params.n = 12;
        %params.Vs = [165.552;106.29;171.479217603912;136.298969072165;139.823;293.267;240.069;342.924648786718;287.885514018692;159.251968503937;399.232262210797;521.991728395062;62326]; %cm3/mol
        params.Vs = [165.550724637681;106.274509803922;171.479217603912;136.298969072165;142.2;285.558638083228;241.834319526627;342.924648786718;287.885514018692;159.251968503937;399.232262210797;492.789393939394;62326]; %cm3/mol
        params.HanSolParam = [14.1,0,0;18,1.4,2;16.2,0,0.6;19.6,2,2.9;20.6,0.8,4.7;16.3,0,0;18,0,0.6;17.7,0,0;17.4,0.1,1.1;19.7,1,2;18.3,0,0;18.7,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [42.00;27.75;1.31;0.3;0.067;0.975;0.1;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [2.04357212833801;0.897528151921969;1.41688074042455;1.21192731038566;0.87463956108188;1.30469280273796;1.27148292183791;1.34918590552474;1.31889268777699;1.24234733094848;1.38608094232405;1.46761444694738];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [0.0382465541337784;0.727735998210257;0.00789547587979925;0.00655969385000662;0.0295234986593422;0.00387261977266118;0.00813801592322982;0.00147867962025794;0.00283721132942447;0.101809993788825;0.00068260813161763;0.00013289924521644]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [0.085362082;0.875671066;0.040986893;0.069035742;0.035283528;0.00760617;0.021697915;0.024962898;0.004882497;0.637351482;0.240391994;280.91]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [0]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.221;0.579;0.321;0;0.748];
            params.bs = [6880;0.973;15800;31.5;255]; %torr^-1 mixed exp and genome
            params.ks = [0.0525;0.00744;0.171;0.538;10.7];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.9500096480641;109.920896713247;17.9512856335608;5.59517317369966;1.36250783890009]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.7676366509148;89.2859269313939;13.4406834054801;4.56821100305585;0.741321469339973]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [89.6901094516630;531.240864816618;62.8682639525013;26.6727574346793;6.35724827979464]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [66.6422526266453;431.513339457563;47.0714159032490;21.7771246968768;3.45888992433533]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.Ch = [0.77;0.68;0.73;0.71;1.01];
            params.bs = [0.59;0.68;5.32;6.55;161.73];
            params.chis_im = [0.73;1.38;0.71;0.57;1.02];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [24.4320861716655;106.715078468667;16.5111963633180;6.11999794342432;1.16129980370385]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.1537212921656;86.6819229437563;12.3624442465211;6.11968813199987;0.631846701207788]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.5760737647037;278.508561986389;7.59319690360824;1.54591628595111;2.06900952589494]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.60134564606960;226.225366152481;5.68526164101519;1.26217215750965;1.12571916140146]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.0661335514905027;0.385089938367672;0.157540194288305;0.218081362599799;0.403204539056517;0.187579254487108;0.197833750525596;0.174871347561278;0.18340051124887;0.207424326859053;0.165145630611683;0.145938384674287];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
    elseif contains(sysInfo.memID,'Ducky-9-12comp')

        params.lmem = 0.3; %thickness of active membrane layer um
        params.compID = struct('IOC', 1, 'TOL', 2, 'BCH', 3, 'TET', 4, 'MNP', 5, 'ICT', 6, 'TPB', 7, 'PRI', 8, 'DDB', 9, 'PYR', 10, 'DCS', 11, 'SQU', 12);
        params.n = 12;
        %params.Vs = [165.552;106.29;171.479217603912;136.298969072165;139.823;293.267;240.069;342.924648786718;287.885514018692;159.251968503937;399.232262210797;521.991728395062;62326]; %cm3/mol
        params.Vs = [165.550724637681;106.274509803922;171.479217603912;136.298969072165;142.2;285.558638083228;241.834319526627;342.924648786718;287.885514018692;159.251968503937;399.232262210797;492.789393939394;62326]; %cm3/mol
        params.HanSolParam = [14.1,0,0;18,1.4,2;16.2,0,0.6;19.6,2,2.9;20.6,0.8,4.7;16.3,0,0;18,0,0.6;17.7,0,0;17.4,0.1,1.1;19.7,1,2;18.3,0,0;18.7,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [42.00;27.75;1.31;0.3;0.067;0.975;0.1;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [2.21045310205603;0.943517165208983;1.5479847136113;1.32054465889681;0.896048860308936;1.40619440453814;1.36127668008067;1.46664183218029;1.4254547101019;1.19256824081476;1.51696697018483;1.6286925665955];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [0.0674849840374727;0.700357039109492;0.0174629897628376;0.00786624151961294;0.0216485385041308;0.00510704710655926;0.0110845198533018;0.00187175471781255;0.0036919816787286;0.0688556104711481;0.000836704742837739;0.000152463677708746]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [0.0516883894180568;0.294694479293164;0.0220878782614502;0.0227191981841919;0.00842655384041737;0.00881259882032173;0.00851106681695843;0.0275290991521341;0.00214109530916375;0.229426193695797;0.302253264378026;546.180260184209]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [0]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.221;0.579;0.321;0;0.748];
            params.bs = [6880;0.973;15800;31.5;255]; %torr^-1 mixed exp and genome
            params.ks = [0.0525;0.00744;0.171;0.538;10.7];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.9500096480641;109.920896713247;17.9512856335608;5.59517317369966;1.36250783890009]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.7676366509148;89.2859269313939;13.4406834054801;4.56821100305585;0.741321469339973]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [89.6901094516630;531.240864816618;62.8682639525013;26.6727574346793;6.35724827979464]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [66.6422526266453;431.513339457563;47.0714159032490;21.7771246968768;3.45888992433533]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.Ch = [0.77;0.68;0.73;0.71;1.01];
            params.bs = [0.59;0.68;5.32;6.55;161.73];
            params.chis_im = [0.73;1.38;0.71;0.57;1.02];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [24.4320861716655;106.715078468667;16.5111963633180;6.11999794342432;1.16129980370385]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.1537212921656;86.6819229437563;12.3624442465211;6.11968813199987;0.631846701207788]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.5760737647037;278.508561986389;7.59319690360824;1.54591628595111;2.06900952589494]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.60134564606960;226.225366152481;5.68526164101519;1.26217215750965;1.12571916140146]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.185517077041459;0.351898890724704;0.129638609821303;0.182922042703839;0.386226787721773;0.160128477546228;0.171606863441495;0.146150645232613;0.155501124896967;0.225231234023961;0.135649951207285;0.115456679958606];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
    elseif contains(sysInfo.memID,'Ducky-10-12comp')

        params.lmem = 0.3; %thickness of active membrane layer um
        params.compID = struct('IOC', 1, 'TOL', 2, 'BCH', 3, 'TET', 4, 'MNP', 5, 'ICT', 6, 'TPB', 7, 'PRI', 8, 'DDB', 9, 'PYR', 10, 'DCS', 11, 'SQU', 12);
        params.n = 12;
        %params.Vs = [165.552;106.29;171.479217603912;136.298969072165;139.823;293.267;240.069;342.924648786718;287.885514018692;159.251968503937;399.232262210797;521.991728395062;62326]; %cm3/mol
        params.Vs = [165.550724637681;106.274509803922;171.479217603912;136.298969072165;142.2;285.558638083228;241.834319526627;342.924648786718;287.885514018692;159.251968503937;399.232262210797;492.789393939394;62326]; %cm3/mol
        params.HanSolParam = [14.1,0,0;18,1.4,2;16.2,0,0.6;19.6,2,2.9;20.6,0.8,4.7;16.3,0,0;18,0,0.6;17.7,0,0;17.4,0.1,1.1;19.7,1,2;18.3,0,0;18.7,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [42.00;27.75;1.31;0.3;0.067;0.975;0.1;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [2.0189611779746;0.88500025855316;1.38886535404286;1.20444563823788;0.864048239556851;1.26671653714554;1.2374073894869;1.30591512203598;1.27923475255981;1.23813055026012;1.33836762032382;1.40994487328983];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [0.0181612945570541;0.586916032655397;0.00502135841113587;0.0062521763386608;0.0318213119773944;0.00384675178935032;0.0076043450308349;0.00158970474376224;0.0028912741223966;0.104678211598224;0.000781892276890587;0.000174043440555237]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [0.027439867;0.44340236;0.01620723;0.049349073;0.026070178;0.004873231;0.01937131;0.021374502;0.004645998;0.563244552;0.252825573;329.92]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [0]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.221;0.579;0.321;0;0.748];
            params.bs = [6880;0.973;15800;31.5;255]; %torr^-1 mixed exp and genome
            params.ks = [0.0525;0.00744;0.171;0.538;10.7];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.9500096480641;109.920896713247;17.9512856335608;5.59517317369966;1.36250783890009]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.7676366509148;89.2859269313939;13.4406834054801;4.56821100305585;0.741321469339973]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [89.6901094516630;531.240864816618;62.8682639525013;26.6727574346793;6.35724827979464]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [66.6422526266453;431.513339457563;47.0714159032490;21.7771246968768;3.45888992433533]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.Ch = [0.77;0.68;0.73;0.71;1.01];
            params.bs = [0.59;0.68;5.32;6.55;161.73];
            params.chis_im = [0.73;1.38;0.71;0.57;1.02];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [24.4320861716655;106.715078468667;16.5111963633180;6.11999794342432;1.16129980370385]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.1537212921656;86.6819229437563;12.3624442465211;6.11968813199987;0.631846701207788]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.5760737647037;278.508561986389;7.59319690360824;1.54591628595111;2.06900952589494]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.60134564606960;226.225366152481;5.68526164101519;1.26217215750965;1.12571916140146]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.0682550985592174;0.394864011954694;0.16443948123369;0.220809453153037;0.411981598474435;0.199363432720581;0.209108835803151;0.187214843210977;0.195377693471756;0.208861149209475;0.177858721147124;0.159214100374249];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
    elseif contains(sysInfo.memID,'Mallard-1-12comp')

        params.lmem = 0.3; %thickness of active membrane layer um
        params.compID = struct('IOC', 1, 'TOL', 2, 'BCH', 3, 'TET', 4, 'MNP', 5, 'ICT', 6, 'TPB', 7, 'PRI', 8, 'DDB', 9, 'PYR', 10, 'DCS', 11, 'SQU', 12);
        params.n = 12;
        %params.Vs = [165.552;106.29;171.479217603912;136.298969072165;139.823;293.267;240.069;342.924648786718;287.885514018692;159.251968503937;399.232262210797;521.991728395062;62326]; %cm3/mol
        params.Vs = [165.550724637681;106.274509803922;171.479217603912;136.298969072165;142.2;285.558638083228;241.834319526627;342.924648786718;287.885514018692;159.251968503937;399.232262210797;492.789393939394;62326]; %cm3/mol
        params.HanSolParam = [14.1,0,0;18,1.4,2;16.2,0,0.6;19.6,2,2.9;20.6,0.8,4.7;16.3,0,0;18,0,0.6;17.7,0,0;17.4,0.1,1.1;19.7,1,2;18.3,0,0;18.7,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [42.00;27.75;1.31;0.3;0.067;0.975;0.1;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [2.38911007669624;1.04817290866032;1.6920661097019;1.41682661110672;0.917369921713943;1.46000289978999;1.4289289432823;1.50128800608237;1.47322057853681;1.15341436944333;1.53524579859453;1.60949213250038];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [27.9930409301753;79.5992309914253;9.39651401870162;1.72748177456251;2.11143305990951;1.57328680014304;3.08426476500285;0.656920662575918;1.18652811136436;1.52037649705238;0.325662873979109;0.0736253523194922]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [6.42221735734438;9.76553408068794;1.54556046470333;0.627189819071333;0.141424938970443;1.09459636102226;0.302404326272558;1.32788052045262;0.0735815615156396;0.941852727968384;7.27510054776659;4164.84650348289]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [0]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.221;0.579;0.321;0;0.748];
            params.bs = [6880;0.973;15800;31.5;255]; %torr^-1 mixed exp and genome
            params.ks = [0.0525;0.00744;0.171;0.538;10.7];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.9500096480641;109.920896713247;17.9512856335608;5.59517317369966;1.36250783890009]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.7676366509148;89.2859269313939;13.4406834054801;4.56821100305585;0.741321469339973]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [89.6901094516630;531.240864816618;62.8682639525013;26.6727574346793;6.35724827979464]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [66.6422526266453;431.513339457563;47.0714159032490;21.7771246968768;3.45888992433533]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.Ch = [0.77;0.68;0.73;0.71;1.01];
            params.bs = [0.59;0.68;5.32;6.55;161.73];
            params.chis_im = [0.73;1.38;0.71;0.57;1.02];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [24.4320861716655;106.715078468667;16.5111963633180;6.11999794342432;1.16129980370385]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.1537212921656;86.6819229437563;12.3624442465211;6.11968813199987;0.631846701207788]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.5760737647037;278.508561986389;7.59319690360824;1.54591628595111;2.06900952589494]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.60134564606960;226.225366152481;5.68526164101519;1.26217215750965;1.12571916140146]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.043084430919867;0.289470519820425;0.10562321273822;0.157553172770336;0.370271302276289;0.147609754041274;0.154684234867734;0.138818738266047;0.144722321945759;0.240641900114799;0.132067112484353;0.118652069707637];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
    elseif contains(sysInfo.memID,'SBAD-1-ML')
      % Lively/Ronita's 9-comp data TOL/MCH/1MN/DCN/NOC/IOC/TBB/TPB/ICE/SBAD-1
        params.lmem = 0.3; % thickness of active membrane layer um
        params.compID = struct('TOL', 1, 'MCH', 2, 'MNP', 3, 'DEC', 4, 'NOC', 5, 'IOC', 6, 'TBB', 7, 'TPB', 8, 'ICT', 9);
        params.n = 9; 
        params.Vs = [106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;62326]; %cm3/mol
        params.HanSolParam = [18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0]; %[delD,delP,delH;...]
        params.psat = [28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458];%torr
        
        if contains(sysInfo.memPhaseModel,'F-H')    
%             params.chis_im = [0.783093904;0.768167694;0.737833686;0.686575528;1.010762328;0.934616296;0.74904465;0.710264759;0.677519155]; %_wo_SBAD
            params.chis_im = [0.999013681;1.649668087;1.281579544;3.74659782;1.68870698;2.232767509;2.37490061;3.005702275;3.859089665]; %_partially_w_SBAD
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [0]; %um2/s FH Exp Based tweaked for single comp MS diffs 
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %Mixed Exp (LOW ERROR) Based FH Tweaked to match single comp FH um2/s   
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0
                    params.diffs_im = [186.6955732;321.1679179;10.50535274;29593390.95;384825.6833;298.1088;10.34303218;1.228239816;5.485]; %um2/s FICKS FH-BC tweaked for single comp 20bar
                elseif sysInfo.lowDiffErrorBar == 1    
                    params.diffs_im = [0]; %um2/s FICKS (LOW ERROR) FH-BC tweaked for single comp 20 bar
                end
            end
        end
    %  FFV Diff Model
        params.unitActPhis = [0.420;0.1;0.55;0.0248;0.1569;0.02;0.0429;0.0288;0.0153]; %used with FFV calcs. Get these values from unit activity exp sorption volume fractions 
        params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1]; 
    end
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%custom parameter sets based on inputs
    compID = [];
    for i = 1:sysInfo.n
        compID_i = getfield(params.compID,sysInfo.mixID(i));
        compID = [compID;compID_i];
    end
    if contains(sysInfo.memPhaseModel,'F-H')||contains(sysInfo.memPhaseModel,'FH-LM')
        params.chis = [0.1*ones(sysInfo.n),params.chis_im(compID,end);params.chis_im([compID;1],end).']; %note chi_ji = chi_ij
    end
    params.diffs = [100000*ones(sysInfo.n),params.diffs_im(compID,end)]; %um2/s  
    params.Vs = [params.Vs(compID);params.Vs(end)]; %cm3/mol
    if sysInfo.crossDiffFudge == 1
        for k = 1:length(sysInfo.crossDiffSpecs)/2
            j = sysInfo.crossDiffSpecs(2*k);
            i = sysInfo.crossDiffSpecs(2*k-1);
            params.diffs(i,j) = sysInfo.crossDiffVals(k);
            params.diffs(j,i) = params.diffs(i,j)*params.Vs(j)/params.Vs(i);
        end
    end
    if contains(sysInfo.memPhaseModel,'DSM')||contains(sysInfo.memPhaseModel,'FH-LM')
        params.bs = params.bs(compID);
        params.Ch = params.Ch(compID);
    end
    if contains(sysInfo.memPhaseModel,'DSM')
        params.ks = params.ks(compID);
    end
    params.HanSolParam = params.HanSolParam(compID,:);
    params.psat = params.psat(compID);
    params.Bffv = params.Bffv(compID);
    params.unitActPhis = params.unitActPhis(compID);
%------------------------------------------------------------------------------------------------------------------------------------% 

%------------------------------------------------------------------------------------------------------------------------------------% 
%carry over needed sysInfo to params
    params.T = sysInfo.T;
    params.R = sysInfo.R;
    params.mixID = sysInfo.mixID;
    if contains(sysInfo.memID,'SBAD1')
        params.memID = 1;
    elseif contains(sysInfo.memID,'PIM1')
        params.memID = 2;
    end
    params.compID = cell2struct(mat2cell(1:sysInfo.n,1,ones(sysInfo.n,1)),params.mixID,2);
    params.Pu = sysInfo.Pu;
    params.Pd = sysInfo.Pd;
    params.yf = sysInfo.yf;
    params.n = sysInfo.n;
    if contains(sysInfo.memPhaseModel,'F-H')
        params.memPhaseModel = 1;
    elseif contains(sysInfo.memPhaseModel,'DSM')
        params.memPhaseModel = 2;
    elseif contains(sysInfo.memPhaseModel,'FH-LM')
        params.memPhaseModel = 3;
    end
    if contains(sysInfo.diffModel,'NoCoupling')
        params.diffModel = 1;
    elseif contains(sysInfo.diffModel,'Vignes')
        params.diffModel = 2;
    elseif contains(sysInfo.diffModel,'Darken')
        params.diffModel = 3;
    elseif contains(sysInfo.diffModel,'Fudge')
        params.diffModel = 4;
    end
    if contains(sysInfo.swlDiffModel,'none')
        params.swlDiffModel = 1;
    elseif contains(sysInfo.swlDiffModel,'FFV-UAV')
        params.swlDiffModel = 2;
    elseif contains(sysInfo.swlDiffModel,'Avg-Diff')
        params.swlDiffModel = 3;
    end
    if contains(sysInfo.numMethod,'FullDis')
        params.numNodes = sysInfo.numNodes; 
    elseif contains(sysInfo.numMethod,'MultShootAlg')
        params.numShootPoints = sysInfo.numShootPoints; 
        params.casADi = sysInfo.casADi;
    end
    params.solverSpec = sysInfo.solverSpec;
    params.iterDetail = sysInfo.iterDetail;
    params.crossDiffFudge = sysInfo.crossDiffFudge;
    params.memPhaseModel_FicksOG = sysInfo.memPhaseModel_FicksOG;
    params.noThermoCoupling = sysInfo.noThermoCoupling;
    if sysInfo.crossDiffFudge == 1
        params.crossDiffSpecs = sysInfo.crossDiffSpecs;
        params.crossDiffVals = sysInfo.crossDiffVals;
    end
    params.pervapMode = sysInfo.pervapMode;
    params.currentStateLit_eqnSetup = sysInfo.currentStateLit_eqnSetup;
    params.noGammaFugacityODEs = sysInfo.noGammaFugacityODEs;
    params.nodalGuessApprox = sysInfo.nodalGuessApprox;
%------------------------------------------------------------------------------------------------------------------------------------% 

end

