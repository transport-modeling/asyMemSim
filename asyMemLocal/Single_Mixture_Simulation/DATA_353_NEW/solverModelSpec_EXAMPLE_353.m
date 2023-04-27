function [sysInfo] = solverModelSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%sorption/diffusion model and thermodynamic assumpstions spec         
  %sorption membrane phase model
    sysInfo.memPhaseModel = "F-H"; %other options: "FH-LM", "F-H", and "DSM"
    
  %specify thermodynamic assumptions
    sysInfo.noThermoCoupling = 0;
    sysInfo.memPhaseModel_FicksOG = 0;
    
  %specify diffusional relationships    
    %sysInfo.diffModel = "Vignes";  %Vignes, NoCoupling
    %sysInfo.swlDiffModel = "none"; % none, Avg-Diff, FFV-UAV
    
  %specify diffusional relationships    
%     sysInfo.diffModel = "Vignes";  %Vignes, NoCoupling
%     sysInfo.swlDiffModel = "FFV-UAV"; % none, Avg-Diff, FFV-UAV
%     sysInfo.diffModel = "NoCoupling";  %Vignes, NoCoupling
%     sysInfo.swlDiffModel = "none"; % none, Avg-Diff, FFV-UAV
   
  %specify diffusional relationships    
    sysInfo.diffModel = "NoCoupling";  %Vignes, NoCoupling
    sysInfo.swlDiffModel = "Avg-Diff"; % none, Avg-Diff, FFV-UAV
%------------------------------------------------------------------------------------------------------------------------------------%

%------------------------------------------------------------------------------------------------------------------------------------%
%%numerical method sepcs
%   sysInfo.numMethod = "NormShootAlg"; %classical shooting algorithm (SA) [Recommended] % "NormShootAlg" % "MultShootAlg"
%  %solver specs
%    sysInfo.solverSpec = "levenberg-marquardt"; %default, use for FUD
%    sysInfo.iterDetail = 1; %0 = all main solver output will be suppressed 
%  %initial guess specs
%    sysInfo.initGuessApprox = 0; %0 = use feed mol fraction vector and small total flux value (use if Forward Euler IVP approx fails)
%  %sys of eqns spec
%    sysInfo.currentStateLit_eqnSetup = 0; %default == 0 (no difference to genreal simulations)
%    sysInfo.noGammaFugacityODEs = 0; % ONLY USE FOR FULL MS and SA or FUD method (FICK == 0 and noThermo == 0)! option keep MS system interms of fugacity gradients and do not make chain rule trick to get gamma matrix
%------------------------------------------------------------------------------------------------------------------------------------%    
%%%%numerical method sepcs
   sysInfo.numMethod = "NormShootAlg"; %classical shooting algorithm (SA) [Recommended]
%      sysInfo.numMethod = "FullDis"; 
%      sysInfo.numNodes = 30;
%    sysInfo.numShootPoints = 5;
%    sysInfo.casADi = 0;
%  %solver specs
    sysInfo.solverSpec = "trust-region-dogleg"; %default, use for FUD
    sysInfo.iterDetail = 1; %0 = all main solver output will be suppressed 
  %initial guess specs
    sysInfo.initGuessApprox = 0; %0 = use feed mol fraction vector and small total flux value (use if Forward Euler IVP approx fails)
    sysInfo.nodalGuessApprox = 0; %only applicable to FUD methods N>1
%     sysInfo.customInitGuessApprox = []; %use this if don't want custom guess
%     sysInfo.customInitGuessApprox = [1E-6*ones(354,1)];

    %FH, no therm/diff coupling solution 353 comp
    sysInfo.customInitGuessApprox = [0.00745228683667722;0.00710447574522774;0.0100793787909057;0.0200732534733599;0.0150985228695196;0.0228285921339456;0.0231295973998379;0.0245869674519531;0.0233453120534681;0.0300948994850963;0.00207202398326659;0.0276165026510413;0.0157806477348252;0.00341609268004130;0.0215595137025277;0.0331273285159087;0.00684454845186293;0.0305785954175874;0.0255628108462630;0.00775004688290199;0.0269403222989641;0.0106584640756630;0.0214301193281699;0.00116368386656430;0.0255026306753704;0.00921786717359600;0.0250761891531512;0.0131771084055386;0.0116605998956014;0.00684297731418701;0.0181662693552335;0.0105249445908403;0.0226824402156632;0.0150151742162061;0.00451799138818464;0.0127466982174138;0.0110160910309282;0.0115820353271740;0.0201344061836589;0.00686785067720369;0.00606616620454329;0.00109230491487016;0.0213235593261281;0.0123978338515162;0.00507475273498767;0.0121005384282479;0.00861982514508099;0.00600883403344950;0.000649258856751091;0.00245654176957172;0.0131713346296499;0.0129919611221530;0.00674334160304698;0.00749599522228951;0.00466368384762315;0.00192698822753527;0.00949443034755164;0.000521743140584855;0.0104972119884634;0.00676212307288674;0.00714130252994344;0.00731918504364453;0.00388212810683900;0.000798982885256687;0.00231693737239739;0.0100634005624401;0.00246172660797283;0.00653498834414013;0.0135981691483719;0.0103747574452850;0.000159663338294977;0.000234597219676554;0.00112069700537534;0.00105723239377570;0.0103489981280019;0.00216242871322943;0.00339301269191308;0.00540217404413937;0.00661852459927102;0.00511888824087733;2.58734842601458e-05;0.00173597173609066;0.000164667906260771;0.000489374273783960;0.00909087754146999;1.37036324107372e-06;0.00167506149267605;0.00105578669495503;0.00477934067850205;0.00637282919490835;0.00275549779608594;9.70898600922976e-05;0.00230390857868533;0.000188566753750127;0.00128526459366772;0.000907373917782946;0.00590206691273826;0.00317737805440875;0.00390755576444233;1.22514422892831e-05;0.00357796228908139;8.73738623406837e-05;0.000876998353581156;0.000220875745517521;0.00160306576372206;0.000839444995099591;0.00351704241706447;0.00182960068817332;0.000917165276207333;1.45626612636664e-05;4.90917179889385e-05;0.000183346736274291;0.00273752464461490;0.00214679144143635;0.00131275682778654;0.000647172733488292;0.000545200298206960;0.000853701548891125;0.00235136176458551;0.00104666990516575;1.31074274915693e-05;0.000823258746602476;0.00105510586346645;0.00184444155394751;0.000446095604117033;0.000283219472949784;0.000224786598240703;0.000801351774497087;0.00105700234939088;0.00476374942188283;0.000108113695740048;6.67843990842956e-05;3.34093362712485e-06;9.39366890273948e-05;4.50926591169820e-05;6.40199768706352e-05;0.00182668567032213;0.000298426807872472;0.00104089662610996;0.000257131476105438;0.000155573550396407;8.72877761145183e-05;9.19791729338316e-06;0.000201887536223741;1.56611585404451e-05;3.84953591721305e-05;5.08095696694494e-05;7.61047098708063e-05;1.61318429142404e-06;2.61489090214988e-06;0.000716596612638214;0.00228837090428397;0.000146387203886076;0.000588497038515923;0.000615743183546057;9.35298147276111e-05;1.73545519318348e-05;3.20309353087132e-05;0.000122233885324990;0.000191864622698001;9.12645024599356e-05;1.75595344157361e-05;0.000217245497285108;0.000776215788542522;6.55886058305868e-05;1.67928993071686e-07;0.000234388555770506;0.000190812058522146;0.000230563402512803;0.000109608855048917;1.85212889055778e-05;7.91455754302700e-05;3.83421340475896e-05;0.00124339108508453;0.00154634817057911;3.15759848184428e-06;9.86585535337794e-06;7.03385049612071e-07;0.000123385345117885;0.000243331135424834;0.000102575498885457;1.51603440643189e-05;0.000258497022397348;3.06226140788178e-05;0.00207236891861512;0.000848084374743889;8.05953169083668e-05;3.57163775140938e-06;3.70211153590982e-07;6.63510488199081e-05;0.000140935434572200;8.61906424519293e-05;1.08226707246076e-05;0.000121843498490581;5.57183363790613e-05;8.88215751472585e-05;1.86860786978360e-07;0.000446996706740822;4.91822947970988e-05;8.47186246244094e-06;4.32615264259043e-06;1.89038062760471e-07;3.47470321551629e-05;1.93575828351250e-05;7.38991761927615e-05;6.77764317515167e-05;7.14209283342681e-06;4.55479156830761e-05;1.13787985169789e-05;7.38759969130821e-09;1.55071106781256e-05;0.000196625248555010;1.33822222329965e-07;6.41140314251167e-05;7.80637520362413e-05;6.46486415921552e-05;2.77676724731740e-06;3.15039850250178e-05;0.000223298373316959;9.95371777787343e-06;5.12063890460446e-05;4.57109539316775e-06;6.93067007529269e-06;8.02382000882014e-05;4.89047329488472e-09;1.48511713232252e-07;7.82456746258151e-08;5.02357175680272e-05;8.57726031449632e-05;9.94290956107116e-05;2.57236606511393e-06;1.65581586274251e-05;1.11948766002597e-05;0.000185661646706254;1.30941551959113e-05;0.000135284723246861;1.35752653232933e-05;8.65073990067365e-07;4.37810499475298e-05;8.59132489436212e-07;4.42622118486430e-05;1.50266835928105e-09;0.000260012867597891;0.000114158975730884;0.000236351025894634;3.67170795709129e-05;4.26927280453637e-07;2.60505421735332e-05;6.68615181291135e-05;1.41645098778404e-06;1.99737619014956e-05;7.65309354710336e-07;5.94526010557143e-06;2.77606272536293e-08;5.87861520689208e-05;9.41937073343677e-05;0.00120767627211885;4.65380762047948e-06;3.61411779416714e-06;1.23272139721287e-05;1.15945183203988e-06;1.56891663608031e-06;3.57922697676472e-05;1.55839940244143e-07;2.37909058646013e-08;7.48832093867137e-05;0.000111322470409357;0.000250466974820217;1.45358928515167e-06;4.55760865775274e-06;7.83218021598163e-05;3.30338269319084e-06;8.19121317810079e-07;3.73678257334186e-07;4.28659459204176e-07;2.02291470410858e-05;2.44596582816189e-07;1.64326741614290e-08;0.000138807050251996;1.47714526639121e-05;1.05713457500254e-05;7.92242292979902e-06;4.30623120060820e-07;3.86967129205532e-06;1.09894078840092e-05;2.56268020578617e-07;1.15422394571707e-08;0.000110276897113469;9.11663590928740e-05;3.40432508504123e-06;7.96992174321052e-07;1.99154362652304e-06;5.59564416930820e-06;2.27113673880725e-07;8.33766070830388e-09;4.52877977730820e-06;3.43306462670493e-06;8.51736094746276e-08;2.30774105404911e-07;2.56213380631473e-06;4.94552220550088e-05;4.16659061802438e-06;7.26049157409827e-05;2.35098427472375e-07;1.24065042965521e-07;9.50299280795821e-07;1.75330810235687e-09;9.07901093951117e-07;1.15801432515131e-08;1.19274544119631e-06;0.000107461089238167;8.45960431643857e-08;6.49536184555127e-08;1.62755579730880e-07;1.84310790724681e-09;3.29873177606614e-06;3.58955414348863e-08;1.58943257645004e-06;1.44468152597272e-07;3.34271842902535e-06;1.20319997219380e-08;1.07240619262116e-06;8.98225463427665e-08;5.03795530348365e-08;7.21710598988521e-07;5.54280899280090e-08;2.32330755646837e-05;2.52225102279478e-09;6.09148316848037e-08;1.10904336717556e-08;3.82447930871310e-08;3.55747490444160e-09;1.06641693518445e-08;4.35459913047382e-08;1.11632906317441e-05;4.05115375698199e-10;2.02487011427277e-08;1.26201636274175e-08;1.44221718416015e-10;7.78156108826184e-09;6.47415101100170e-08;2.49125300913441e-07;5.76683272866572e-09;2.12718288860245e-07;9.86771339586729e-11;1.57829313716549e-09;9.64734904694425e-10;1.27628501331965e-09;1.95238956893060e-10;1.45045576210716e-09;1.97764431716984e-10;1.38584239201156e-10;5.35012556762882e-10;0.00653551950408940]; %set to match your experimental compostion number plus total flux as [n+1 by 1] vector for custom initial guess
    %FH, thermo coupling, no diff coupling, AVG diff  
    %sysInfo.customInitGuessApprox = [0.00700021465839948;0.00668330444787366;0.00953529756943288;0.0185098105586110;0.0138876470830337;0.0210575776027983;0.0215269095731036;0.0237037337852809;0.0209892497006164;0.0271778049510618;0.00189106428722857;0.0253113930381079;0.0139511114662134;0.00324349254517706;0.0189187556672386;0.0292387076463496;0.00616141406950540;0.0265089882159825;0.0218834321319690;0.00726560352295155;0.0231946067082811;0.00931323318674140;0.0191402639919608;0.000996960141726465;0.0212813171750745;0.00847825264731507;0.0210291694111094;0.0113967266661704;0.00992983431800927;0.00602272110786529;0.0148221737447971;0.00947415914249580;0.0185457116252720;0.0125794138663891;0.00384216374231523;0.0111808375047776;0.00879542192644333;0.0101920009725337;0.0160788730189044;0.00619251032137718;0.00519009556501873;0.000934580572396866;0.0167094756706043;0.0106435934451062;0.00409124829017475;0.00946997745741944;0.00748945738028162;0.00563004214292137;0.000525408203825204;0.00216117587617330;0.0101536271937600;0.0108619069007202;0.00561934318312422;0.00577332257496314;0.00434221149083650;0.00173074218926443;0.00817844378800167;0.000419530852928905;0.00800623782196523;0.00624169846475361;0.00580733801926007;0.00653382487976068;0.00296718850461143;0.000723765060056830;0.00222704890164152;0.00855010312953253;0.00199330426757381;0.00500369484522997;0.0107595922921838;0.00911955499495840;0.000182550954714310;0.000241628005096114;0.000865357255278305;0.00112157079363608;0.00877737783168208;0.00305072045900827;0.00268016265567126;0.00471245563033596;0.00510757084837330;0.00670803473654592;2.20688795495897e-05;0.00205288251587025;0.000132977631307712;0.000637964487377366;0.00792583278778177;0.000231102790013277;0.00352112663816819;0.000901308405620335;0.00467468923477662;0.00481846121543735;0.00407233852477451;8.02873571725955e-05;0.00174081406614770;0.000249666160186803;0.00115218676334971;0.00167518633890946;0.00544387155772693;0.00315657235433058;0.00456751597352877;0.000785931082741583;0.00618493956497891;7.00273114459956e-05;0.00132117927741070;0.000166395720855877;0.00175181502122360;0.00239893567092107;0.00355535741738600;0.00230678142732001;0.00448284784606874;0.00120735903865056;0.000729302064290545;0.000322896470406314;0.00213075805749913;0.00163307564544938;0.00196530079357899;0.00288372057556134;0.000925801647599293;0.000980887252239368;0.00216624344939394;0.00213823475076568;0.00148498369254314;0.00173863997449904;0.000823320325586913;0.00140457817746195;0.00101852455581305;0.00314073120830956;0.00153954386329101;0.00187129232788435;0.00143114545788606;0.00428842344820736;0.000310271578159210;0.00247440653700531;0.000531520526234322;0.000243172141108648;0.00140237987265216;5.28047877587876e-05;0.00137357252790604;0.00110339828416907;0.000925946136724451;0.00275756008328615;0.000261653859250778;0.000277897142026744;8.54447729875388e-05;0.00378557172160377;0.000142502826101894;0.00192010348254693;4.49856996495831e-05;0.000245477850166388;7.62873813123469e-05;0.000595571014659349;0.000657774931530877;0.00170858770287304;0.000861166169139747;0.000553512597803703;0.00259721979144585;0.000228946794237304;0.000236531691853627;7.07283943112936e-05;0.00310865609431583;0.000779178341457155;7.96768453256846e-05;0.00128331622194754;0.000167842396441244;0.000578671407750193;7.19859557170718e-05;5.84205170146311e-05;0.00205284242563170;0.00197897584303158;0.00123563857859271;0.000353041174627308;0.000383595602081229;0.00265726171472008;0.000197369764100185;0.00107460730384875;0.00123170285780717;0.000341055617630533;1.40459703846262e-05;0.000387309446520718;0.00149466503054625;0.00159476946720275;0.000439594561463457;0.000490757441709613;0.00170207886833787;0.000279774433627973;0.00181073164106601;0.000716068083868677;0.000159713020170433;0.000552911845760243;0.000325897279635577;0.00104850497229922;0.00110954991621964;0.000497211679735011;0.000566277212405567;0.00104431828718133;0.000659974119495132;0.000776583244725449;0.000221089624722226;0.000405857713154663;4.56031307933918e-05;2.40503642260036e-05;0.000969920465583635;0.000270959083212758;0.000692722789406161;3.05494445859148e-05;0.000703217723915862;0.000532697253879825;0.000617010215625314;0.000522417973541227;0.000967182376408378;4.13492341246995e-05;0.000180702868714223;0.000192685745865001;0.000230485474662994;0.000256723754701907;8.22107962185138e-05;0.000240344075389793;0.000929566279464185;0.000778973212129066;0.000450544729292127;0.000118141415311161;0.000551346491376343;0.000648280925596687;0.000110281755375479;0.00128198819567569;3.48646266896022e-05;5.21063447436114e-05;0.000186753867066811;0.000241509239538850;0.000105425333805809;0.000109295289908947;0.00125139863251794;0.000503319083894840;0.000175005235237858;0.000497649680448714;0.000193556974699940;0.000469796061304478;0.000311309410834953;0.000199499299382627;0.000986709857604294;0.000410076906992789;0.000282337583102986;4.73610174908724e-06;0.000381288647926084;0.000591054853304931;0.000304791312607235;0.000812264915492885;1.60778066217077e-05;0.000537483994520790;0.000246729149888458;0.000520686790832443;0.000651114670036863;0.000491190666975122;5.10653105150388e-05;0.000112077021412385;0.000435024847470370;0.000166289676051615;0.00191978958350512;0.000156973364184117;0.000174328280584665;0.000366872956939815;0.000664176026507346;0.000397237857475985;0.000411164481172434;0.000132602642938281;0.000122414697204507;0.000189914498938070;0.000236837696993887;0.00221798482631082;7.92171497484239e-05;0.000295447245822207;0.000155009774640809;0.000145915630565090;0.000337204081307095;0.000324637385914556;0.000168738325129678;0.000309511246934151;0.000273292811889381;0.000107943759730218;0.000354945805963687;4.58234906438634e-05;0.000720140413041645;0.000523501722780725;0.000254763286911353;0.000448876718547862;0.000225389208824876;0.000378717889228869;9.59631407340235e-05;0.000337690078968610;0.000347131044074485;0.000323626841393982;0.000676138110456555;0.000369027331611996;0.000155787627363174;0.000455825069203793;8.59067919578052e-05;0.000598257680289768;0.000643024795793112;0.000103960544164189;0.000315813527016886;9.83567791624510e-05;0.000300894234799480;0.000761857670226732;0.000419318290795613;0.000416710064739657;0.000255932130520646;5.11260121827610e-05;1.46573862417726e-05;0.000323921222125104;0.000172705610557170;0.000301092071877981;0.000757362850931956;0.000221381642697412;0.000205296359824801;1.24492083478496e-05;1.94714212610913e-05;0.000254646983829628;1.93322709852647e-05;0.000550919254327745;0.000569596062421865;0.000363316851642807;4.82850552605130e-05;0.000507168840231435;0.000542435919090495;0.000316203812862132;0.000467239446929837;0.000515760917484162;0.000317192609456145;2.42829440661298e-05;0.000865784098788154;0.000156672522925042;0.000811522516817815;4.93620241786174e-05;0.000508782273633690;0.000297078111644117;0.000694009091521223;2.73541456644541e-05;1.30341639064926e-05;0.000584234921782940;5.78094845793004e-06;0.00101366869944246;0.000429715630527698;0.000479304568553797;0.000265899812246312;0.000109642445195247;1.95153674984605e-05;0.000163555409235084;0.000184867687009093;0.000536591904657493;6.65230386800104e-05;0.000389226639277500;5.77741928853502e-05;4.40105146774973e-05;0.000348032444287264;0.0518752066856201];
   %sys of eqns spec
   sysInfo.currentStateLit_eqnSetup = 0; %default == 0 (no difference to genreal simulations)
    sysInfo.noGammaFugacityODEs = 0; % ONLY USE FOR FULL MS and SA or FUD method (FICK == 0 and noThermo == 0)! option keep MS system interms of fugacity gradients and do not make chain rule trick to get gamma matrix
%------------------------------------------------------------------------------------------------------------------------------------% 