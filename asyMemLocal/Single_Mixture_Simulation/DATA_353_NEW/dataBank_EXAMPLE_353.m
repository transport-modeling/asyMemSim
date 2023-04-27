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
            params.chis_im = [2.088385991;0.997650735;1.270355078;1.585798617;1.054722721;2.429960601;2.123303626;2.027403056;1.540047314;1.388287198;1.746780602;2.253015493];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [0.710648;1.96393;0.70721;0.724879;1.3407;0.151848;0.207784;0.0875696;0.136232;0.764001;0.0641025;0.017023]; %um2/s FH EXP (tweaked @20 bar)
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
    params.unitActPhis = [0.062463317;0.317583193;0.198194358;0.122744578;0.286070751;0.041021076;0.05976592;0.067518609;0.131145279;0.164585795;0.097933971;0.050851089];
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
            params.chis_im = [2.018348657;1.070831164;1.31820903;1.551404898;1.185693209;2.290521959;2.145593377;1.915831949;1.547014129;1.568547616;1.590057936;2.102244064];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [0.479496;1.96137;0.297021;0.349838;0.417755;0.0268073;0.0383033;0.0129004;0.0197363;0.13106;0.00877726;0.0016039]; %um2/s FH EXP (tweaked @20 bar)
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
    params.unitActPhis = [0.068308898;0.277930569;0.183598985;0.128995928;0.227842905;0.048563115;0.058114095;0.078051713;0.129821704;0.125832725;0.12199637;0.061376336];
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
            params.chis_im = [1.980433833;1.004749641;1.233865051;1.581165104;1.094328502;2.283907264;2.107020162;1.922460713;1.503007723;1.441550102;1.639222217;2.149373942];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [3.65569;7.71184;4.61961;4.77021;9.17513;1.49512;2.10336;0.888779;1.5091;8.50179;0.660814;0.183536]; %um2/s FH EXP (tweaked @20 bar)
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
    params.unitActPhis = [0.071737092;0.313420118;0.210327551;0.123564874;0.266594166;0.048958002;0.061006776;0.07737465;0.138466762;0.151761163;0.113748579;0.057839157];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
    
    
    elseif contains(sysInfo.memID,'Ducky-10-12comp')

        params.lmem = 0.3; %thickness of active membrane layer um
        params.compID = struct("1", 1, "2", 2, "3", 3, "4", 4, "5", 5, "6", 6, "7", 7, "8", 8, "9", 9, "10", 10, "11", 11, "12", 12);
        params.n = 12;
        %params.Vs = [165.552;106.29;171.479217603912;136.298969072165;139.823;293.267;240.069;342.924648786718;287.885514018692;159.251968503937;399.232262210797;521.991728395062;62326]; %cm3/mol
        params.Vs = [165.550724637681;106.274509803922;171.479217603912;136.298969072165;142.2;285.558638083228;241.834319526627;342.924648786718;287.885514018692;159.251968503937;399.232262210797;492.789393939394;62326]; %cm3/mol
        params.HanSolParam = [14.1,0,0;18,1.4,2;16.2,0,0.6;19.6,2,2.9;20.6,0.8,4.7;16.3,0,0;18,0,0.6;17.7,0,0;17.4,0.1,1.1;19.7,1,2;18.3,0,0;18.7,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [42.00;27.75;1.31;0.3;0.067;0.975;0.1;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [1.980433833;1.004749641;1.233865051;1.581165104;1.094328502;2.283907264;2.107020162;1.922460713;1.503007723;1.441550102;1.639222217;2.149373942];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [3.65569;7.71184;4.61961;4.77021;9.17513;1.49512;2.10336;0.888779;1.5091;8.50179;0.660814;0.183536]; %um2/s FH EXP (tweaked @20 bar)
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
    params.unitActPhis = [0.071737092;0.313420118;0.210327551;0.123564874;0.266594166;0.048958002;0.061006776;0.07737465;0.138466762;0.151761163;0.113748579;0.057839157];
    params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1];
    
  elseif contains(sysInfo.memID,'CrudeOil-DUCKY9')

        params.lmem = .3; %nm %0.3; %thickness of active membrane layer um
        %params.compID = struct('SOL1', 1, 'SOL2', 2, 'SOL3', 3, 'SOL4', 4, 'SOL5', 5, 'SOL6', 6, 'SOL7', 7, 'SOL8', 8, 'SOL9', 9, 'SOL10', 10, 'SOL11', 11, 'SOL12', 12, 'SOL13', 13, 'SOL14', 14, 'SOL15', 15, 'SOL16', 16, 'SOL17', 17, 'SOL18', 18, 'SOL19', 19, 'SOL20', 20, 'SOL21', 21, 'SOL22', 22, 'SOL23', 23, 'SOL24', 24, 'SOL25', 25, 'SOL26', 26, 'SOL27', 27, 'SOL28', 28, 'SOL29', 29, 'SOL30', 30, 'SOL31', 31, 'SOL32', 32, 'SOL33', 33, 'SOL34', 34, 'SOL35', 35, 'SOL36', 36, 'SOL37', 37, 'SOL38', 38, 'SOL39', 39, 'SOL40', 40, 'SOL41', 41, 'SOL42', 42, 'SOL43', 43, 'SOL44', 44, 'SOL45', 45, 'SOL46', 46, 'SOL47', 47, 'SOL48', 48, 'SOL49', 49, 'SOL50', 50, 'SOL51', 51, 'SOL52', 52, 'SOL53', 53, 'SOL54', 54, 'SOL55', 55, 'SOL56', 56, 'SOL57', 57, 'SOL58', 58, 'SOL59', 59, 'SOL60', 60, 'SOL61', 61, 'SOL62', 62, 'SOL63', 63, 'SOL64', 64, 'SOL65', 65, 'SOL66', 66, 'SOL67', 67, 'SOL68', 68, 'SOL69', 69, 'SOL70', 70, 'SOL71', 71, 'SOL72', 72, 'SOL73', 73, 'SOL74', 74, 'SOL75', 75, 'SOL76', 76, 'SOL77', 77, 'SOL78', 78, 'SOL79', 79, 'SOL80', 80, 'SOL81', 81, 'SOL82', 82, 'SOL83', 83, 'SOL84', 84, 'SOL85', 85, 'SOL86', 86, 'SOL87', 87, 'SOL88', 88, 'SOL89', 89, 'SOL90', 90, 'SOL91', 91, 'SOL92', 92, 'SOL93', 93, 'SOL94', 94, 'SOL95', 95, 'SOL96', 96, 'SOL97', 97, 'SOL98', 98, 'SOL99', 99, 'SOL100', 100, 'SOL101', 101, 'SOL102', 102, 'SOL103', 103, 'SOL104', 104, 'SOL105', 105, 'SOL106', 106, 'SOL107', 107, 'SOL108', 108, 'SOL109', 109, 'SOL110', 110, 'SOL111', 111, 'SOL112', 112, 'SOL113', 113, 'SOL114', 114, 'SOL115', 115, 'SOL116', 116, 'SOL117', 117, 'SOL118', 118, 'SOL119', 119, 'SOL120', 120, 'SOL121', 121, 'SOL122', 122, 'SOL123', 123, 'SOL124', 124, 'SOL125', 125, 'SOL126', 126, 'SOL127', 127, 'SOL128', 128, 'SOL129', 129, 'SOL130', 130, 'SOL131', 131, 'SOL132', 132, 'SOL133', 133, 'SOL134', 134, 'SOL135', 135, 'SOL136', 136, 'SOL137', 137, 'SOL138', 138, 'SOL139', 139, 'SOL140', 140, 'SOL141', 141, 'SOL142', 142, 'SOL143', 143, 'SOL144', 144, 'SOL145', 145, 'SOL146', 146, 'SOL147', 147, 'SOL148', 148, 'SOL149', 149, 'SOL150', 150, 'SOL151', 151, 'SOL152', 152, 'SOL153', 153, 'SOL154', 154, 'SOL155', 155, 'SOL156', 156, 'SOL157', 157, 'SOL158', 158, 'SOL159', 159, 'SOL160', 160, 'SOL161', 161, 'SOL162', 162, 'SOL163', 163, 'SOL164', 164, 'SOL165', 165, 'SOL166', 166, 'SOL167', 167, 'SOL168', 168, 'SOL169', 169, 'SOL170', 170, 'SOL171', 171, 'SOL172', 172, 'SOL173', 173, 'SOL174', 174, 'SOL175', 175, 'SOL176', 176, 'SOL177', 177, 'SOL178', 178, 'SOL179', 179, 'SOL180', 180, 'SOL181', 181, 'SOL182', 182, 'SOL183', 183, 'SOL184', 184, 'SOL185', 185, 'SOL186', 186, 'SOL187', 187, 'SOL188', 188, 'SOL189', 189, 'SOL190', 190, 'SOL191', 191, 'SOL192', 192, 'SOL193', 193, 'SOL194', 194, 'SOL195', 195, 'SOL196', 196, 'SOL197', 197, 'SOL198', 198, 'SOL199', 199, 'SOL200', 200, 'SOL201', 201, 'SOL202', 202, 'SOL203', 203, 'SOL204', 204, 'SOL205', 205, 'SOL206', 206, 'SOL207', 207, 'SOL208', 208, 'SOL209', 209, 'SOL210', 210, 'SOL211', 211, 'SOL212', 212, 'SOL213', 213, 'SOL214', 214, 'SOL215', 215, 'SOL216', 216, 'SOL217', 217, 'SOL218', 218, 'SOL219', 219, 'SOL220', 220, 'SOL221', 221, 'SOL222', 222, 'SOL223', 223, 'SOL224', 224, 'SOL225', 225, 'SOL226', 226, 'SOL227', 227, 'SOL228', 228, 'SOL229', 229, 'SOL230', 230, 'SOL231', 231, 'SOL232', 232, 'SOL233', 233, 'SOL234', 234, 'SOL235', 235, 'SOL236', 236, 'SOL237', 237, 'SOL238', 238, 'SOL239', 239, 'SOL240', 240, 'SOL241', 241, 'SOL242', 242, 'SOL243', 243, 'SOL244', 244, 'SOL245', 245, 'SOL246', 246, 'SOL247', 247, 'SOL248', 248, 'SOL249', 249, 'SOL250', 250, 'SOL251', 251, 'SOL252', 252, 'SOL253', 253, 'SOL254', 254, 'SOL255', 255, 'SOL256', 256, 'SOL257', 257, 'SOL258', 258, 'SOL259', 259, 'SOL260', 260, 'SOL261', 261, 'SOL262', 262, 'SOL263', 263, 'SOL264', 264, 'SOL265', 265, 'SOL266', 266, 'SOL267', 267, 'SOL268', 268, 'SOL269', 269, 'SOL270', 270, 'SOL271', 271, 'SOL272', 272, 'SOL273', 273, 'SOL274', 274, 'SOL275', 275, 'SOL276', 276, 'SOL277', 277, 'SOL278', 278, 'SOL279', 279, 'SOL280', 280, 'SOL281', 281, 'SOL282', 282, 'SOL283', 283, 'SOL284', 284, 'SOL285', 285, 'SOL286', 286, 'SOL287', 287, 'SOL288', 288, 'SOL289', 289, 'SOL290', 290, 'SOL291', 291, 'SOL292', 292, 'SOL293', 293, 'SOL294', 294, 'SOL295', 295, 'SOL296', 296, 'SOL297', 297, 'SOL298', 298, 'SOL299', 299, 'SOL300', 300, 'SOL301',301, 'SOL302', 302, 'SOL303', 303, 'SOL304', 304, 'SOL305', 305, 'SOL306', 306, 'SOL307', 307, 'SOL308', 308, 'SOL309', 309, 'SOL310', 310, 'SOL311', 311, 'SOL312', 312, 'SOL313', 313, 'SOL314', 314, 'SOL315', 315, 'SOL316', 316, 'SOL317', 317, 'SOL318', 318, 'SOL319', 319, 'SOL320', 320, 'SOL321', 321, 'SOL322', 322, 'SOL323', 323, 'SOL324', 324, 'SOL325', 325, 'SOL326', 326, 'SOL327', 327, 'SOL328', 328, 'SOL329', 329, 'SOL330', 330, 'SOL331', 331, 'SOL332', 332, 'SOL333', 333, 'SOL334', 334, 'SOL335', 335, 'SOL336', 336, 'SOL337', 337, 'SOL338', 338, 'SOL339', 339, 'SOL340', 340, 'SOL341', 341, 'SOL342', 342, 'SOL343', 343, 'SOL344', 344, 'SOL345', 345, 'SOL346', 346, 'SOL347', 347, 'SOL348', 348, 'SOL349', 349, 'SOL350', 350, 'SOL351', 351, 'SOL352', 352, 'SOL353', 353, 'SOL354', 354, 'SOL355', 355, 'SOL356', 356, 'SOL357', 357, 'SOL358', 358, 'SOL359', 359, 'SOL360', 360, 'SOL361', 361, 'SOL362', 362, 'SOL363', 363, 'SOL364', 364, 'SOL365', 365, 'SOL366', 366, 'SOL367', 367, 'SOL368', 368, 'SOL369', 369, 'SOL370', 370, 'SOL371', 371);
        params.compID = struct('SOL1', 1, 'SOL2', 2, 'SOL3', 3, 'SOL4', 4, 'SOL5', 5, 'SOL6', 6, 'SOL7', 7, 'SOL8', 8, 'SOL9', 9, 'SOL10', 10, 'SOL11', 11, 'SOL12', 12, 'SOL13', 13, 'SOL14', 14, 'SOL15', 15, 'SOL16', 16, 'SOL17', 17, 'SOL18', 18, 'SOL19', 19, 'SOL20', 20, 'SOL21', 21, 'SOL22', 22, 'SOL23', 23, 'SOL24', 24, 'SOL25', 25, 'SOL26', 26, 'SOL27', 27, 'SOL28', 28, 'SOL29', 29, 'SOL30', 30, 'SOL31', 31, 'SOL32', 32, 'SOL33', 33, 'SOL34', 34, 'SOL35', 35, 'SOL36', 36, 'SOL37', 37, 'SOL38', 38, 'SOL39', 39, 'SOL40', 40, 'SOL41', 41, 'SOL42', 42, 'SOL43', 43, 'SOL44', 44, 'SOL45', 45, 'SOL46', 46, 'SOL47', 47, 'SOL48', 48, 'SOL49', 49, 'SOL50', 50, 'SOL51', 51, 'SOL52', 52, 'SOL53', 53, 'SOL54', 54, 'SOL55', 55, 'SOL56', 56, 'SOL57', 57, 'SOL58', 58, 'SOL59', 59, 'SOL60', 60, 'SOL61', 61, 'SOL62', 62, 'SOL63', 63, 'SOL64', 64, 'SOL65', 65, 'SOL66', 66, 'SOL67', 67, 'SOL68', 68, 'SOL69', 69, 'SOL70', 70, 'SOL71', 71, 'SOL72', 72, 'SOL73', 73, 'SOL74', 74, 'SOL75', 75, 'SOL76', 76, 'SOL77', 77, 'SOL78', 78, 'SOL79', 79, 'SOL80', 80, 'SOL81', 81, 'SOL82', 82, 'SOL83', 83, 'SOL84', 84, 'SOL85', 85, 'SOL86', 86, 'SOL87', 87, 'SOL88', 88, 'SOL89', 89, 'SOL90', 90, 'SOL91', 91, 'SOL92', 92, 'SOL93', 93, 'SOL94', 94, 'SOL95', 95, 'SOL96', 96, 'SOL97', 97, 'SOL98', 98, 'SOL99', 99, 'SOL100', 100, 'SOL101', 101, 'SOL102', 102, 'SOL103', 103, 'SOL104', 104, 'SOL105', 105, 'SOL106', 106, 'SOL107', 107, 'SOL108', 108, 'SOL109', 109, 'SOL110', 110, 'SOL111', 111, 'SOL112', 112, 'SOL113', 113, 'SOL114', 114, 'SOL115', 115, 'SOL116', 116, 'SOL117', 117, 'SOL118', 118, 'SOL119', 119, 'SOL120', 120, 'SOL121', 121, 'SOL122', 122, 'SOL123', 123, 'SOL124', 124, 'SOL125', 125, 'SOL126', 126, 'SOL127', 127, 'SOL128', 128, 'SOL129', 129, 'SOL130', 130, 'SOL131', 131, 'SOL132', 132, 'SOL133', 133, 'SOL134', 134, 'SOL135', 135, 'SOL136', 136, 'SOL137', 137, 'SOL138', 138, 'SOL139', 139, 'SOL140', 140, 'SOL141', 141, 'SOL142', 142, 'SOL143', 143, 'SOL144', 144, 'SOL145', 145, 'SOL146', 146, 'SOL147', 147, 'SOL148', 148, 'SOL149', 149, 'SOL150', 150, 'SOL151', 151, 'SOL152', 152, 'SOL153', 153, 'SOL154', 154, 'SOL155', 155, 'SOL156', 156, 'SOL157', 157, 'SOL158', 158, 'SOL159', 159, 'SOL160', 160, 'SOL161', 161, 'SOL162', 162, 'SOL163', 163, 'SOL164', 164, 'SOL165', 165, 'SOL166', 166, 'SOL167', 167, 'SOL168', 168, 'SOL169', 169, 'SOL170', 170, 'SOL171', 171, 'SOL172', 172, 'SOL173', 173, 'SOL174', 174, 'SOL175', 175, 'SOL176', 176, 'SOL177', 177, 'SOL178', 178, 'SOL179', 179, 'SOL180', 180, 'SOL181', 181, 'SOL182', 182, 'SOL183', 183, 'SOL184', 184, 'SOL185', 185, 'SOL186', 186, 'SOL187', 187, 'SOL188', 188, 'SOL189', 189, 'SOL190', 190, 'SOL191', 191, 'SOL192', 192, 'SOL193', 193, 'SOL194', 194, 'SOL195', 195, 'SOL196', 196, 'SOL197', 197, 'SOL198', 198, 'SOL199', 199, 'SOL200', 200, 'SOL201', 201, 'SOL202', 202, 'SOL203', 203, 'SOL204', 204, 'SOL205', 205, 'SOL206', 206, 'SOL207', 207, 'SOL208', 208, 'SOL209', 209, 'SOL210', 210, 'SOL211', 211, 'SOL212', 212, 'SOL213', 213, 'SOL214', 214, 'SOL215', 215, 'SOL216', 216, 'SOL217', 217, 'SOL218', 218, 'SOL219', 219, 'SOL220', 220, 'SOL221', 221, 'SOL222', 222, 'SOL223', 223, 'SOL224', 224, 'SOL225', 225, 'SOL226', 226, 'SOL227', 227, 'SOL228', 228, 'SOL229', 229, 'SOL230', 230, 'SOL231', 231, 'SOL232', 232, 'SOL233', 233, 'SOL234', 234, 'SOL235', 235, 'SOL236', 236, 'SOL237', 237, 'SOL238', 238, 'SOL239', 239, 'SOL240', 240, 'SOL241', 241, 'SOL242', 242, 'SOL243', 243, 'SOL244', 244, 'SOL245', 245, 'SOL246', 246, 'SOL247', 247, 'SOL248', 248, 'SOL249', 249, 'SOL250', 250, 'SOL251', 251, 'SOL252', 252, 'SOL253', 253, 'SOL254', 254, 'SOL255', 255, 'SOL256', 256, 'SOL257', 257, 'SOL258', 258, 'SOL259', 259, 'SOL260', 260, 'SOL261', 261, 'SOL262', 262, 'SOL263', 263, 'SOL264', 264, 'SOL265', 265, 'SOL266', 266, 'SOL267', 267, 'SOL268', 268, 'SOL269', 269, 'SOL270', 270, 'SOL271', 271, 'SOL272', 272, 'SOL273', 273, 'SOL274', 274, 'SOL275', 275, 'SOL276', 276, 'SOL277', 277, 'SOL278', 278, 'SOL279', 279, 'SOL280', 280, 'SOL281', 281, 'SOL282', 282, 'SOL283', 283, 'SOL284', 284, 'SOL285', 285, 'SOL286', 286, 'SOL287', 287, 'SOL288', 288, 'SOL289', 289, 'SOL290', 290, 'SOL291', 291, 'SOL292', 292, 'SOL293', 293, 'SOL294', 294, 'SOL295', 295, 'SOL296', 296, 'SOL297', 297, 'SOL298', 298, 'SOL299', 299, 'SOL300', 300, 'SOL301',301, 'SOL302', 302, 'SOL303', 303, 'SOL304', 304, 'SOL305', 305, 'SOL306', 306, 'SOL307', 307, 'SOL308', 308, 'SOL309', 309, 'SOL310', 310, 'SOL311', 311, 'SOL312', 312, 'SOL313', 313, 'SOL314', 314, 'SOL315', 315, 'SOL316', 316, 'SOL317', 317, 'SOL318', 318, 'SOL319', 319, 'SOL320', 320, 'SOL321', 321, 'SOL322', 322, 'SOL323', 323, 'SOL324', 324, 'SOL325', 325, 'SOL326', 326, 'SOL327', 327, 'SOL328', 328, 'SOL329', 329, 'SOL330', 330, 'SOL331', 331, 'SOL332', 332, 'SOL333', 333, 'SOL334', 334, 'SOL335', 335, 'SOL336', 336, 'SOL337', 337, 'SOL338', 338, 'SOL339', 339, 'SOL340', 340, 'SOL341',341, 'SOL342',342, 'SOL343',343, 'SOL344',344,'SOL345',345, 'SOL346',346, 'SOL347',347, 'SOL348',348, 'SOL349',349, 'SOL350',350, 'SOL351',351, 'SOL352',352, 'SOL353',353);
        params.n = 353;
        params.Vs = xlsread('dataMolarVolume');
        params.HanSolParam = xlsread('dataHansen'); %[delD,delP,delH;...]
        params.psat = xlsread('dataVaporPressure');
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = xlsread('dataChi');
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = xlsread('dataDiffusion');%um2/s FH EXP (tweaked @20 bar)
                %trying nm2/s to maybe fix scaling issuess
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
    params.unitActPhis = xlsread('dataPhi');
    params.Bffv = xlsread('dataFreeVolume');
    
    elseif contains(sysInfo.memID,'9-comp-permPrediction')

        params.lmem = 0.3; %thickness of active membrane layer um
        %params.compID = struct('SOL1', 1, 'SOL2', 2, 'SOL3', 3, 'SOL4', 4, 'SOL5', 5, 'SOL6', 6, 'SOL7', 7, 'SOL8', 8, 'SOL9', 9, 'SOL10', 10, 'SOL11', 11, 'SOL12', 12, 'SOL13', 13, 'SOL14', 14, 'SOL15', 15, 'SOL16', 16, 'SOL17', 17, 'SOL18', 18, 'SOL19', 19, 'SOL20', 20, 'SOL21', 21, 'SOL22', 22, 'SOL23', 23, 'SOL24', 24, 'SOL25', 25, 'SOL26', 26, 'SOL27', 27, 'SOL28', 28, 'SOL29', 29, 'SOL30', 30, 'SOL31', 31, 'SOL32', 32, 'SOL33', 33, 'SOL34', 34, 'SOL35', 35, 'SOL36', 36, 'SOL37', 37, 'SOL38', 38, 'SOL39', 39, 'SOL40', 40, 'SOL41', 41, 'SOL42', 42, 'SOL43', 43, 'SOL44', 44, 'SOL45', 45, 'SOL46', 46, 'SOL47', 47, 'SOL48', 48, 'SOL49', 49, 'SOL50', 50, 'SOL51', 51, 'SOL52', 52, 'SOL53', 53, 'SOL54', 54, 'SOL55', 55, 'SOL56', 56, 'SOL57', 57, 'SOL58', 58, 'SOL59', 59, 'SOL60', 60, 'SOL61', 61, 'SOL62', 62, 'SOL63', 63, 'SOL64', 64, 'SOL65', 65, 'SOL66', 66, 'SOL67', 67, 'SOL68', 68, 'SOL69', 69, 'SOL70', 70, 'SOL71', 71, 'SOL72', 72, 'SOL73', 73, 'SOL74', 74, 'SOL75', 75, 'SOL76', 76, 'SOL77', 77, 'SOL78', 78, 'SOL79', 79, 'SOL80', 80, 'SOL81', 81, 'SOL82', 82, 'SOL83', 83, 'SOL84', 84, 'SOL85', 85, 'SOL86', 86, 'SOL87', 87, 'SOL88', 88, 'SOL89', 89, 'SOL90', 90, 'SOL91', 91, 'SOL92', 92, 'SOL93', 93, 'SOL94', 94, 'SOL95', 95, 'SOL96', 96, 'SOL97', 97, 'SOL98', 98, 'SOL99', 99, 'SOL100', 100, 'SOL101', 101, 'SOL102', 102, 'SOL103', 103, 'SOL104', 104, 'SOL105', 105, 'SOL106', 106, 'SOL107', 107, 'SOL108', 108, 'SOL109', 109, 'SOL110', 110, 'SOL111', 111, 'SOL112', 112, 'SOL113', 113, 'SOL114', 114, 'SOL115', 115, 'SOL116', 116, 'SOL117', 117, 'SOL118', 118, 'SOL119', 119, 'SOL120', 120, 'SOL121', 121, 'SOL122', 122, 'SOL123', 123, 'SOL124', 124, 'SOL125', 125, 'SOL126', 126, 'SOL127', 127, 'SOL128', 128, 'SOL129', 129, 'SOL130', 130, 'SOL131', 131, 'SOL132', 132, 'SOL133', 133, 'SOL134', 134, 'SOL135', 135, 'SOL136', 136, 'SOL137', 137, 'SOL138', 138, 'SOL139', 139, 'SOL140', 140, 'SOL141', 141, 'SOL142', 142, 'SOL143', 143, 'SOL144', 144, 'SOL145', 145, 'SOL146', 146, 'SOL147', 147, 'SOL148', 148, 'SOL149', 149, 'SOL150', 150, 'SOL151', 151, 'SOL152', 152, 'SOL153', 153, 'SOL154', 154, 'SOL155', 155, 'SOL156', 156, 'SOL157', 157, 'SOL158', 158, 'SOL159', 159, 'SOL160', 160, 'SOL161', 161, 'SOL162', 162, 'SOL163', 163, 'SOL164', 164, 'SOL165', 165, 'SOL166', 166, 'SOL167', 167, 'SOL168', 168, 'SOL169', 169, 'SOL170', 170, 'SOL171', 171, 'SOL172', 172, 'SOL173', 173, 'SOL174', 174, 'SOL175', 175, 'SOL176', 176, 'SOL177', 177, 'SOL178', 178, 'SOL179', 179, 'SOL180', 180, 'SOL181', 181, 'SOL182', 182, 'SOL183', 183, 'SOL184', 184, 'SOL185', 185, 'SOL186', 186, 'SOL187', 187, 'SOL188', 188, 'SOL189', 189, 'SOL190', 190, 'SOL191', 191, 'SOL192', 192, 'SOL193', 193, 'SOL194', 194, 'SOL195', 195, 'SOL196', 196, 'SOL197', 197, 'SOL198', 198, 'SOL199', 199, 'SOL200', 200, 'SOL201', 201, 'SOL202', 202, 'SOL203', 203, 'SOL204', 204, 'SOL205', 205, 'SOL206', 206, 'SOL207', 207, 'SOL208', 208, 'SOL209', 209, 'SOL210', 210, 'SOL211', 211, 'SOL212', 212, 'SOL213', 213, 'SOL214', 214, 'SOL215', 215, 'SOL216', 216, 'SOL217', 217, 'SOL218', 218, 'SOL219', 219, 'SOL220', 220, 'SOL221', 221, 'SOL222', 222, 'SOL223', 223, 'SOL224', 224, 'SOL225', 225, 'SOL226', 226, 'SOL227', 227, 'SOL228', 228, 'SOL229', 229, 'SOL230', 230, 'SOL231', 231, 'SOL232', 232, 'SOL233', 233, 'SOL234', 234, 'SOL235', 235, 'SOL236', 236, 'SOL237', 237, 'SOL238', 238, 'SOL239', 239, 'SOL240', 240, 'SOL241', 241, 'SOL242', 242, 'SOL243', 243, 'SOL244', 244, 'SOL245', 245, 'SOL246', 246, 'SOL247', 247, 'SOL248', 248, 'SOL249', 249, 'SOL250', 250, 'SOL251', 251, 'SOL252', 252, 'SOL253', 253, 'SOL254', 254, 'SOL255', 255, 'SOL256', 256, 'SOL257', 257, 'SOL258', 258, 'SOL259', 259, 'SOL260', 260, 'SOL261', 261, 'SOL262', 262, 'SOL263', 263, 'SOL264', 264, 'SOL265', 265, 'SOL266', 266, 'SOL267', 267, 'SOL268', 268, 'SOL269', 269, 'SOL270', 270, 'SOL271', 271, 'SOL272', 272, 'SOL273', 273, 'SOL274', 274, 'SOL275', 275, 'SOL276', 276, 'SOL277', 277, 'SOL278', 278, 'SOL279', 279, 'SOL280', 280, 'SOL281', 281, 'SOL282', 282, 'SOL283', 283, 'SOL284', 284, 'SOL285', 285, 'SOL286', 286, 'SOL287', 287, 'SOL288', 288, 'SOL289', 289, 'SOL290', 290, 'SOL291', 291, 'SOL292', 292, 'SOL293', 293, 'SOL294', 294, 'SOL295', 295, 'SOL296', 296, 'SOL297', 297, 'SOL298', 298, 'SOL299', 299, 'SOL300', 300, 'SOL301',301, 'SOL302', 302, 'SOL303', 303, 'SOL304', 304, 'SOL305', 305, 'SOL306', 306, 'SOL307', 307, 'SOL308', 308, 'SOL309', 309, 'SOL310', 310, 'SOL311', 311, 'SOL312', 312, 'SOL313', 313, 'SOL314', 314, 'SOL315', 315, 'SOL316', 316, 'SOL317', 317, 'SOL318', 318, 'SOL319', 319, 'SOL320', 320, 'SOL321', 321, 'SOL322', 322, 'SOL323', 323, 'SOL324', 324, 'SOL325', 325, 'SOL326', 326, 'SOL327', 327, 'SOL328', 328, 'SOL329', 329, 'SOL330', 330, 'SOL331', 331, 'SOL332', 332, 'SOL333', 333, 'SOL334', 334, 'SOL335', 335, 'SOL336', 336, 'SOL337', 337, 'SOL338', 338, 'SOL339', 339, 'SOL340', 340, 'SOL341', 341, 'SOL342', 342, 'SOL343', 343, 'SOL344', 344, 'SOL345', 345, 'SOL346', 346, 'SOL347', 347, 'SOL348', 348, 'SOL349', 349, 'SOL350', 350, 'SOL351', 351, 'SOL352', 352, 'SOL353', 353, 'SOL354', 354, 'SOL355', 355, 'SOL356', 356, 'SOL357', 357, 'SOL358', 358, 'SOL359', 359, 'SOL360', 360, 'SOL361', 361, 'SOL362', 362, 'SOL363', 363, 'SOL364', 364, 'SOL365', 365, 'SOL366', 366, 'SOL367', 367, 'SOL368', 368, 'SOL369', 369, 'SOL370', 370, 'SOL371', 371);
        %params.compID = struct('SOL1', 1, 'SOL2', 2, 'SOL3', 3, 'SOL4', 4, 'SOL5', 5, 'SOL6', 6, 'SOL7', 7, 'SOL8', 8, 'SOL9', 9, 'SOL10', 10, 'SOL11', 11, 'SOL12', 12, 'SOL13', 13, 'SOL14', 14, 'SOL15', 15, 'SOL16', 16, 'SOL17', 17, 'SOL18', 18, 'SOL19', 19, 'SOL20', 20, 'SOL21', 21, 'SOL22', 22, 'SOL23', 23, 'SOL24', 24, 'SOL25', 25, 'SOL26', 26, 'SOL27', 27, 'SOL28', 28, 'SOL29', 29, 'SOL30', 30, 'SOL31', 31, 'SOL32', 32, 'SOL33', 33, 'SOL34', 34, 'SOL35', 35, 'SOL36', 36, 'SOL37', 37, 'SOL38', 38, 'SOL39', 39, 'SOL40', 40, 'SOL41', 41, 'SOL42', 42, 'SOL43', 43, 'SOL44', 44, 'SOL45', 45, 'SOL46', 46, 'SOL47', 47, 'SOL48', 48, 'SOL49', 49, 'SOL50', 50, 'SOL51', 51, 'SOL52', 52, 'SOL53', 53, 'SOL54', 54, 'SOL55', 55, 'SOL56', 56, 'SOL57', 57, 'SOL58', 58, 'SOL59', 59, 'SOL60', 60, 'SOL61', 61, 'SOL62', 62, 'SOL63', 63, 'SOL64', 64, 'SOL65', 65, 'SOL66', 66, 'SOL67', 67, 'SOL68', 68, 'SOL69', 69, 'SOL70', 70, 'SOL71', 71, 'SOL72', 72, 'SOL73', 73, 'SOL74', 74, 'SOL75', 75, 'SOL76', 76, 'SOL77', 77, 'SOL78', 78, 'SOL79', 79, 'SOL80', 80, 'SOL81', 81, 'SOL82', 82, 'SOL83', 83, 'SOL84', 84, 'SOL85', 85, 'SOL86', 86, 'SOL87', 87, 'SOL88', 88, 'SOL89', 89, 'SOL90', 90, 'SOL91', 91, 'SOL92', 92, 'SOL93', 93, 'SOL94', 94, 'SOL95', 95, 'SOL96', 96, 'SOL97', 97, 'SOL98', 98, 'SOL99', 99, 'SOL100', 100, 'SOL101', 101, 'SOL102', 102, 'SOL103', 103, 'SOL104', 104, 'SOL105', 105, 'SOL106', 106, 'SOL107', 107, 'SOL108', 108, 'SOL109', 109, 'SOL110', 110, 'SOL111', 111, 'SOL112', 112, 'SOL113', 113, 'SOL114', 114, 'SOL115', 115, 'SOL116', 116, 'SOL117', 117, 'SOL118', 118, 'SOL119', 119, 'SOL120', 120, 'SOL121', 121, 'SOL122', 122, 'SOL123', 123, 'SOL124', 124, 'SOL125', 125, 'SOL126', 126, 'SOL127', 127, 'SOL128', 128, 'SOL129', 129, 'SOL130', 130, 'SOL131', 131, 'SOL132', 132, 'SOL133', 133, 'SOL134', 134, 'SOL135', 135, 'SOL136', 136, 'SOL137', 137, 'SOL138', 138, 'SOL139', 139, 'SOL140', 140, 'SOL141', 141, 'SOL142', 142, 'SOL143', 143, 'SOL144', 144, 'SOL145', 145, 'SOL146', 146, 'SOL147', 147, 'SOL148', 148, 'SOL149', 149, 'SOL150', 150, 'SOL151', 151, 'SOL152', 152, 'SOL153', 153, 'SOL154', 154, 'SOL155', 155, 'SOL156', 156, 'SOL157', 157, 'SOL158', 158, 'SOL159', 159, 'SOL160', 160, 'SOL161', 161, 'SOL162', 162, 'SOL163', 163, 'SOL164', 164, 'SOL165', 165, 'SOL166', 166, 'SOL167', 167, 'SOL168', 168, 'SOL169', 169, 'SOL170', 170, 'SOL171', 171, 'SOL172', 172, 'SOL173', 173, 'SOL174', 174, 'SOL175', 175, 'SOL176', 176, 'SOL177', 177, 'SOL178', 178, 'SOL179', 179, 'SOL180', 180, 'SOL181', 181, 'SOL182', 182, 'SOL183', 183, 'SOL184', 184, 'SOL185', 185, 'SOL186', 186, 'SOL187', 187, 'SOL188', 188, 'SOL189', 189, 'SOL190', 190, 'SOL191', 191, 'SOL192', 192, 'SOL193', 193, 'SOL194', 194, 'SOL195', 195, 'SOL196', 196, 'SOL197', 197, 'SOL198', 198, 'SOL199', 199, 'SOL200', 200, 'SOL201', 201, 'SOL202', 202, 'SOL203', 203, 'SOL204', 204, 'SOL205', 205, 'SOL206', 206, 'SOL207', 207, 'SOL208', 208, 'SOL209', 209, 'SOL210', 210, 'SOL211', 211, 'SOL212', 212, 'SOL213', 213, 'SOL214', 214, 'SOL215', 215, 'SOL216', 216, 'SOL217', 217, 'SOL218', 218, 'SOL219', 219, 'SOL220', 220, 'SOL221', 221, 'SOL222', 222, 'SOL223', 223, 'SOL224', 224, 'SOL225', 225, 'SOL226', 226, 'SOL227', 227, 'SOL228', 228, 'SOL229', 229, 'SOL230', 230, 'SOL231', 231, 'SOL232', 232, 'SOL233', 233, 'SOL234', 234, 'SOL235', 235, 'SOL236', 236, 'SOL237', 237, 'SOL238', 238, 'SOL239', 239, 'SOL240', 240, 'SOL241', 241, 'SOL242', 242, 'SOL243', 243, 'SOL244', 244, 'SOL245', 245, 'SOL246', 246, 'SOL247', 247, 'SOL248', 248, 'SOL249', 249, 'SOL250', 250, 'SOL251', 251, 'SOL252', 252, 'SOL253', 253, 'SOL254', 254, 'SOL255', 255, 'SOL256', 256, 'SOL257', 257, 'SOL258', 258, 'SOL259', 259, 'SOL260', 260, 'SOL261', 261, 'SOL262', 262, 'SOL263', 263, 'SOL264', 264, 'SOL265', 265, 'SOL266', 266, 'SOL267', 267, 'SOL268', 268, 'SOL269', 269, 'SOL270', 270, 'SOL271', 271, 'SOL272', 272, 'SOL273', 273, 'SOL274', 274, 'SOL275', 275, 'SOL276', 276, 'SOL277', 277, 'SOL278', 278, 'SOL279', 279, 'SOL280', 280, 'SOL281', 281, 'SOL282', 282, 'SOL283', 283, 'SOL284', 284, 'SOL285', 285, 'SOL286', 286, 'SOL287', 287, 'SOL288', 288, 'SOL289', 289, 'SOL290', 290, 'SOL291', 291, 'SOL292', 292, 'SOL293', 293, 'SOL294', 294, 'SOL295', 295, 'SOL296', 296, 'SOL297', 297, 'SOL298', 298, 'SOL299', 299, 'SOL300', 300, 'SOL301',301, 'SOL302', 302, 'SOL303', 303, 'SOL304', 304, 'SOL305', 305, 'SOL306', 306, 'SOL307', 307, 'SOL308', 308, 'SOL309', 309, 'SOL310', 310, 'SOL311', 311, 'SOL312', 312, 'SOL313', 313, 'SOL314', 314, 'SOL315', 315, 'SOL316', 316, 'SOL317', 317, 'SOL318', 318, 'SOL319', 319, 'SOL320', 320, 'SOL321', 321, 'SOL322', 322, 'SOL323', 323, 'SOL324', 324, 'SOL325', 325, 'SOL326', 326, 'SOL327', 327, 'SOL328', 328, 'SOL329', 329, 'SOL330', 330, 'SOL331', 331, 'SOL332', 332, 'SOL333', 333, 'SOL334', 334, 'SOL335', 335, 'SOL336', 336, 'SOL337', 337, 'SOL338', 338, 'SOL339', 339, 'SOL340', 340, 'SOL341', 341);
        params.compID = struct('TOL', 1, 'MCH', 2, 'MNP', 3, 'DEC', 4, 'NOC', 5, 'IOC', 6, 'TBB', 7, 'TPB', 8, 'ICT', 9);
        params.n = 9; 
        params.Vs = [106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;62326]; %cm3/mol
        params.HanSolParam = [18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0]; %[delD,delP,delH;...]
        params.psat = [28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458];%torr
        
        if contains(sysInfo.memPhaseModel,'F-H')    
            params.chis_im = [1.036743561;1.117508517;1.15635622;1.677494792;1.179862622;1.559752747;1.348993412;1.351627644;1.803316143];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [10.8545096661386;3.1485072493575;0.0225256250977678;0.121369186763117;11.4744017554721;4.82192410669967;0.535061632422813;0.0230522957116404;0.191054581442974]; %um2/s FH Exp Based tweaked for single comp MS diffs 
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [0]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)   
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
        params.unitActPhis = [0.295532083;0.255997298;0.239437467;0.107791825;0.230088505;0.127443582;0.174923949;0.174205835;0.09068049]; %used with FFV calcs. Get these values from unit activity exp sorption volume fractions 
        params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1]; 
    
     elseif contains(sysInfo.memID,'xyleneOctaneSeparation')

        params.lmem = 0.6; %thickness of active membrane layer um
        params.compID = struct('MXY', 1, 'IOC', 2);
        params.n = 2;
        params.Vs = [123.44;165.55]; %cm3/mol
        params.HanSolParam = [18.4,2.6,2.3;14.1,0,0]; %[delD,delP,delH;...]
        %params.psat = [49.087;28.998;5.25;0.18;0.059;0.0458;0.0352;0.0044;0.000051;0.0000045;0.00000128;0.000000031425]; %torr
        params.psat = [7.675;45.207]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [0.806595362;2.637885821]; % HERE !!!!!
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [14.06462568;7.156630828]; %um2/s FH EXP (tweaked @20 bar)% HERE !!!!!
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
    params.unitActPhis = [0.464480046;0.032089924]; % HERE !!!!!
    params.Bffv = 0.03*[1;1]; 
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

