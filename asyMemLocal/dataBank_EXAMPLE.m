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
            params.chis_im = [0.871;1.672;0.705;2.783;1.163;3.049;1.648;2.5;3.13];
%             params.chis_im = [0.827;1.652;0.678;2.667;1.129;2.949;1.648;2.339;2.912]; % high error
%             params.chis_im = [0.93;1.693;0.743;2.86;1.202;3.241;1.648;2.765;3.61]; %low error
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [16.7701308511911;5.57162654063454;0.0534797474304072;1.66224496710637;5.71465197660387;6.997847822971986;1.40048712268577;0.0691441589002619;0.590373456746892]; %um2/s FH Exp Based tweaked for single comp MS diffs 
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [15.3534759619839;2.86110552131555;0.0425857248048064;1.29285719652783;4.43333987446973;4.217246039008540;0.974251911434123;0.0691441589002619;0.537106678694842]; %Mixed Exp (LOW ERROR) Based FH Tweaked to match single comp FH um2/s   
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0
                    params.diffs_im = [4.21117147213822;3.32944378753008;0.00962069984261326;1.39186144125104;2.36281698371032;6.114523405394954;0.833944291340556;0.0554924921379044;0.521517243958073]; %um2/s FICKS FH-BC tweaked for single comp 20bar
                elseif sysInfo.lowDiffErrorBar == 1    
                    params.diffs_im = [3.85543324275734;1.70971437749573;0.00766092765208338;1.08255889875082;1.83303739099595;3.684911456231122;0.580135159196378;0.0554924921379044;0.474463056537618]; %um2/s FICKS (LOW ERROR) FH-BC tweaked for single comp 20 bar
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.314;0.0730;0.405;0.00759;0.051;0.012;0.0158;0.0137;0.00801];
            params.bs = [0.141;0.18;0;0;0.000019;0.202;0;66.6;67.6]; %torr^-1 mixed exp and genome
            params.ks = [0.0138;0.00121;18.5;0.0263;0.0148;0.000158;0.0525;0.792;0.258]; 
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [13.1543664344202;4.96586274686747;0.0359445909605920;1.57460604635303;6.25103288975875;6.554144601694238;1.34113217627728;0.0592985536611731;0.496597565284541]; %Mixed Exp Based DSM Tweaked to match single comp backup1.22 i-OCT for 20 bar fit but 2.852 for fit of 30bar-60bar %um2/s
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [12.0431528311342;2.55003762676260;0.0286225446480865;1.22469359128711;4.84945601052820;3.949848733459900;0.932961513929973;0.0592985536611731;0.451791769914806]; %Mixed Exp (LOW ERROR) Based DSM Tweaked to match single comp DSM um2/s
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0
                    params.diffs_im = [31.7299700061699;10.5650590499813;0.0775364306435485;1.61776520129334;7.72145931069900;13.749703299546810;1.50019568413218;0.0748186662605057;0.672347124862815]; %um2/s FICKS DSM-BC tweaked for single comp
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [29.0495844110888;5.42530059321879;0.0617419725494922;1.25826182287103;5.99019040975325;8.286245034818150;1.04361438895868;0.0748186662605057;0.611684226373950]; %um2/s FICKS (LOW ERROR) DSM-BC tweaked for single comp
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.chis_im = [1.09;2.32;0.775;3.06;1.25;3.75;1.73;2.86;3.64];
            params.Ch = [0.314;0.0730;0.405;0.00759;0.051;0.0112;0.0158;0.0137;0.00801];
            params.bs = [0.472;0.278;88.8;2.12;0.0164;0.168;0.258;150;110]; %torr^-1 
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [12.6041403148214;5.29540109824826;0.0346308573803066;1.78119158312004;6.37149334526816;6.865419510134641;1.48846000039839;0.0693261706366819;0.596922405281312]; %Mixed Exp Based FHLM Tweaked to match single comp FH-DSM um2/s
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [11.5394070003806;2.71926002341543;0.0275764234694993;1.38537123093526;4.94290739029466;4.137438247819231;1.03545043505932;0.0693261706366819;0.543064744642501];
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [9.87224267781197;9.58130395856575;0.0139037027832745;1.84865074337946;2.98090972577414;12.6311629901020;0.960732596086683;0.0846224274648881;0.866056273551871]; %um2/s FICKS FHLM-BC tweaked for single comp
                elseif sysInfo.lowDiffErrorBar == 1 
                    params.diffs_im = [9.03828610458803;4.92012906077157;0.0110714670307768;1.43783946695811;2.31254431651160;7.612157828467595;0.668335719018764;0.0846224274648881;0.787915857890620]; %um2/s FICKS (LOW ERROR) FHLM-BC tweaked for single comp
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
            params.chis_im = [0.648;0.826;0.642;0.56;0.698];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [50.3055962175403;184.278801966931;35.0494452080733;16.2074097830652;2.80941841324399]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [37.3784609416795;149.684947449161;26.2426049134344;13.2326320429522;1.52856528754996]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [6.83432114205229;45.4102531988667;4.87536536909399;1.31261700013715;0.613893311027835]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [5.07809120717587;36.8855847283164;3.65033701457585;1.07169362670079;0.334010769242320]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.770000000000000;0.679000000000000;0.858000000000000;0.844000000000000;0.812000000000000];
            params.bs = [0.0566;0.595;0.618;0;146]; %torr^-1 mixed exp and genome
            params.ks = [0.0443;0.00401;0.117;0.538;14.3];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.7463328971589;119.665317456880;18.8873864072697;5.60540456372360;1.16983300023063]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.6162991413797;97.2010701363798;14.1415710409058;4.57656448006091;0.636489782934200]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [94.3514434949619;876.702632402083;80.8834710983979;26.6727574346793;5.07969909319987]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [70.1057538162763;712.123076861904;60.5599593245533;21.7771246968768;2.76379327021854]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.chis_im = [0.726;1.65;0.724;0.570;0.891];
            params.Ch = [0.726;0.679;0.858;0.844;0.812];
            params.bs = [0.59;0.726;2.48;2.99;331];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [25.1354668777666;119.057524186868;15.5169176448061;6.84912381815477;1.17373717202742]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.6763528791265;96.7073752417546;11.6179967137487;5.59200615546814;0.638613988202282]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.3078190076714;500.446835859304;8.22030227020068;1.52221964656985;1.39474387010788]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.40202488040216;406.500134202625;6.15479484679679;1.24282490131386;0.758860643305972]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.65;0.46;0.67;0.85;0.6];
    params.Bffv = 0.03*[1;1;1;1;1]; 
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
                params.diffs_im = [926.9722;0;0;0;0;0;0;0;0]; %um2/s FH Exp Based tweaked for single comp MS diffs 
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
    elseif contains(sysInfo.memID,'SBAD-100-C')
      % Lively/Ronita's 9-comp data TOL/MCH/1MN/DCN/NOC/IOC/TBB/TPB/ICE/SBAD-1
        params.lmem = 0.3; % thickness of active membrane layer um
        params.compID = struct('TOL', 1, 'MCH', 2, 'MNP', 3, 'DEC', 4, 'NOC', 5, 'IOC', 6, 'TBB', 7, 'TPB', 8, 'ICT', 9,...
            'TOD', 10, 'MCD', 11, 'MND', 12, 'DED', 13, 'NOD', 14, 'IOD', 15, 'TBD', 16, 'TPD', 17, 'ICD', 18,...
            'TOE', 19, 'MCE', 20, 'MNE', 21, 'DEE', 22, 'NOE', 23, 'IOE', 24, 'TBE', 25, 'TPE', 26, 'ICE', 27);
        params.n = 27; 
        params.Vs = [106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;...
            106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;...
            106.521;128.123;139.823;156.962;163.42;165.552;155.529;240.069;293.267;62326]; %cm3/mol
        params.HanSolParam = [18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0;...
            18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0;...
            18,1.4,2;16,0,1;20.6,0.8,4.7;18,0,0;15.5,0,0;14.1,0,0;17.4,0.1,1.1;18,0,0.6;16.3,0,0]; %[delD,delP,delH;...]
        params.psat = [28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458;...
            28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458;...
            28.998;46.596;0.059;0.975;14.805;49.087;2.115;0.0352;0.0458];%torr
        
        if contains(sysInfo.memPhaseModel,'F-H')    
            params.chis_im = [0.839;1.672;0.705;2.785;1.164;3.046;1.648;2.5;3.128;...
                0.839;1.672;0.705;2.785;1.164;3.046;1.648;2.5;3.128;...
                0.839;1.672;0.705;2.785;1.164;3.046;1.648;2.5;3.128];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [15.4638846409154;5.59021311138405;0.0535253755307862;1.66732636951057;5.74769726115217;6.9787;1.40290524341143;0.0692133718931809;0.589420725298558;...
                    15.4638846409154;5.59021311138405;0.0535253755307862;1.66732636951057;5.74769726115217;6.9787;1.40290524341143;0.0692133718931809;0.589420725298558;...
                    15.4638846409154;5.59021311138405;0.0535253755307862;1.66732636951057;5.74769726115217;6.9787;1.40290524341143;0.0692133718931809;0.589420725298558]; %um2/s FH Exp Based tweaked for single comp MS diffs 
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
            params.chis_im = [1.01;2.32;0.79;3.06;1.25;3.76;1.69;2.66;3.63;...
                1.01;2.32;0.79;3.06;1.25;3.76;1.69;2.66;3.63;...
                1.01;2.32;0.79;3.06;1.25;3.76;1.69;2.66;3.63];
            params.Ch = [0.31;0.07302;0.481;0.00761;0.051;0.0112;0.00821;0.00689;0.00801;...
                0.31;0.07302;0.481;0.00761;0.051;0.0112;0.00821;0.00689;0.00801;...
                0.31;0.07302;0.481;0.00761;0.051;0.0112;0.00821;0.00689;0.00801];
            params.bs = [0.431;0.278;86.057;2.145;0.016;0.157;0.083;509.449;110.020;...
                0.431;0.278;86.057;2.145;0.016;0.157;0.083;509.449;110.020;...
                0.431;0.278;86.057;2.145;0.016;0.157;0.083;509.449;110.020]; %torr^-1 
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [11.4609361017534;5.30042942577299;0.0328169584688351;1.77959539744715;6.38183951490904;6.9251;1.47399503331522;0.0678674014003992;0.593101878582401;...
                    11.4609361017534;5.30042942577299;0.0328169584688351;1.77959539744715;6.38183951490904;6.9251;1.47399503331522;0.0678674014003992;0.593101878582401;...
                    11.4609361017534;5.30042942577299;0.0328169584688351;1.77959539744715;6.38183951490904;6.9251;1.47399503331522;0.0678674014003992;0.593101878582401]; %Mixed Exp Based FHLM Tweaked to match single comp FH-DSM um2/s
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
        params.unitActPhis = [0.420;0.1;0.55;0.0248;0.1569;0.02;0.0429;0.0288;0.0153;...
            0.420;0.1;0.55;0.0248;0.1569;0.02;0.0429;0.0288;0.0153;...
            0.420;0.1;0.55;0.0248;0.1569;0.02;0.0429;0.0288;0.0153]; %used with FFV calcs. Get these values from unit activity exp sorption volume fractions 
        params.Bffv = 0.03*[1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1;1]; 
    elseif contains(sysInfo.memID,'MATRIMID')

      % Lively/Ronita's Exp Data for 5-comp TOL/HEP/PXY/OXY/ICE/PIM-1
        params.lmem = 1.1; %thickness of active membrane layer um
        params.compID = struct('TOL', 1, 'MES', 2, 'TPB', 3);
        params.n = 3;
        params.Vs = [106.521;139.1088;241.8343;71429]; %cm3/mol
        params.HanSolParam = [18,1.4,2;17.8,1,3.1;18,0,0.6]; %[delD,delP,delH;...]
        params.psat = [28.998;2;0.0352]; %torr
        
        if contains(sysInfo.memPhaseModel,'F-H')
            params.chis_im = [0.767;1.445;1.436];
            if sysInfo.lowDiffErrorBar == 0     
                params.diffs_im = [4.22726765798626;1.19271717861478;0.0802916770676725]; %um2/s FH EXP (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [4.22726765798626;1.19271717861478;0.0802916770676725]; %um2/s FH EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [6.83432114205229;45.4102531988667;4.87536536909399]; %um2/s FH-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [5.07809120717587;36.8855847283164;3.65033701457585]; %um2/s FH-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'DSM')
            params.Ch = [0.770000000000000;0.679000000000000;0.858000000000000;0.844000000000000;0.812000000000000];
            params.bs = [0.0566;0.595;0.618;0;146]; %torr^-1 mixed exp and genome
            params.ks = [0.0443;0.00401;0.117;0.538;14.3];
            if sysInfo.lowDiffErrorBar == 0 %noraml diffusivities
                params.diffs_im = [27.7463328971589;119.665317456880;18.8873864072697;5.60540456372360;1.16983300023063]; %um2/s DSM EXP  (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1 %low error diffusivities
                params.diffs_im = [20.6162991413797;97.2010701363798;14.1415710409058;4.57656448006091;0.636489782934200]; %um2/s DSM EXP (LOW ERROR)(tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [94.3514434949619;876.702632402083;80.8834710983979;26.6727574346793;5.07969909319987]; %um2/s DSM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [70.1057538162763;712.123076861904;60.5599593245533;21.7771246968768;2.76379327021854]; %um2/s DSM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        elseif contains(sysInfo.memPhaseModel,'FH-LM')
            params.chis_im = [0.726;1.65;0.724;0.570;0.891];
            params.Ch = [0.726;0.679;0.858;0.844;0.812];
            params.bs = [0.59;0.726;2.48;2.99;331];
            if sysInfo.lowDiffErrorBar == 0
                params.diffs_im = [25.1354668777666;119.057524186868;15.5169176448061;6.84912381815477;1.17373717202742]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            elseif sysInfo.lowDiffErrorBar == 1
                params.diffs_im = [18.6763528791265;96.7073752417546;11.6179967137487;5.59200615546814;0.638613988202282]; %um2/s FH_DSM EXP (LOW ERROR) (tweaked @20 bar)
            end
            if sysInfo.memPhaseModel_FicksOG == 1
                if sysInfo.lowDiffErrorBar == 0 
                    params.diffs_im = [11.3078190076714;500.446835859304;8.22030227020068;1.52221964656985;1.39474387010788]; %um2/s FHLM-BC FICKS  (tweaked @20 bar)
                elseif sysInfo.lowDiffErrorBar == 1
                    params.diffs_im = [8.40202488040216;406.500134202625;6.15479484679679;1.24282490131386;0.758860643305972]; %um2/s FHLM-BC FICKS (LOW ERROR)(tweaked @20 bar)
                end
            end
        end
  %FFV Diff Model
    params.unitActPhis = [0.65;0.46;0.67;0.85;0.6];
    params.Bffv = 0.03*[1;1;1;1;1];         
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
    elseif contains(sysInfo.swlDiffModel,'FFV')
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

