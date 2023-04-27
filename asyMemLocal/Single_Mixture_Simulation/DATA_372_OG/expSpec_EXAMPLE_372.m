function [sysInfo] = expSpec_EXAMPLE()

%------------------------------------------------------------------------------------------------------------------------------------%
%specify single simulation
 
%mixture components, compositions, system parameters
  %SBAD1 Data
    sysInfo.memID = "CrudeOil-SBAD-1";  % polymer spec
    %12 COMP M1
    %sysInfo.mixID = ["1","2","3","4","5","6","7","8","9","10","11","12" ];  % mixture spec
    %sysInfo.mixID = ["IOC","TOL","BCH","TET","MNP","ICT","TPB","PRI","DDB","PYR","DCS","SQU" ];  % mixture spec
    %sysInfo.mixID = ["SOL1","SOL2"];  % mixture spec
    sysInfo.mixID = ["SOL1","SOL2","SOL3","SOL4","SOL5","SOL6","SOL7","SOL8","SOL9","SOL10","SOL11","SOL12","SOL13","SOL14","SOL15","SOL16","SOL17","SOL18","SOL19","SOL20","SOL21","SOL22","SOL23","SOL24","SOL25","SOL26","SOL27","SOL28","SOL29","SOL30","SOL31","SOL32","SOL33","SOL34","SOL35","SOL36","SOL37","SOL38","SOL39","SOL40","SOL41","SOL42","SOL43","SOL44","SOL45","SOL46","SOL47","SOL48","SOL49","SOL50","SOL51","SOL52","SOL53","SOL54","SOL55","SOL56","SOL57","SOL58","SOL59","SOL60","SOL61","SOL62","SOL63","SOL64","SOL65","SOL66","SOL67","SOL68","SOL69","SOL70","SOL71","SOL72","SOL73","SOL74","SOL75","SOL76","SOL77","SOL78","SOL79","SOL80","SOL81","SOL82","SOL83","SOL84","SOL85","SOL86","SOL87","SOL88","SOL89","SOL90","SOL91","SOL92","SOL93","SOL94","SOL95","SOL96","SOL97","SOL98","SOL99","SOL100","SOL101","SOL102","SOL103","SOL104","SOL105","SOL106","SOL107","SOL108","SOL109","SOL110","SOL111","SOL112","SOL113","SOL114","SOL115","SOL116","SOL117","SOL118","SOL119","SOL120","SOL121","SOL122","SOL123","SOL124","SOL125","SOL126","SOL127","SOL128","SOL129","SOL130","SOL131","SOL132","SOL133","SOL134","SOL135","SOL136","SOL137","SOL138","SOL139","SOL140","SOL141","SOL142","SOL143","SOL144","SOL145","SOL146","SOL147","SOL148","SOL149","SOL150","SOL151","SOL152","SOL153","SOL154","SOL155","SOL156","SOL157","SOL158","SOL159","SOL160","SOL161","SOL162","SOL163","SOL164","SOL165","SOL166","SOL167","SOL168","SOL169","SOL170","SOL171","SOL172","SOL173","SOL174","SOL175","SOL176","SOL177","SOL178","SOL179","SOL180","SOL181","SOL182","SOL183","SOL184","SOL185","SOL186","SOL187","SOL188","SOL189","SOL190","SOL191","SOL192","SOL193","SOL194","SOL195","SOL196","SOL197","SOL198","SOL199","SOL200","SOL201","SOL202","SOL203","SOL204","SOL205","SOL206","SOL207","SOL208","SOL209","SOL210","SOL211","SOL212","SOL213","SOL214","SOL215","SOL216","SOL217","SOL218","SOL219","SOL220","SOL221","SOL222","SOL223","SOL224","SOL225","SOL226","SOL227","SOL228","SOL229","SOL230","SOL231","SOL232","SOL233","SOL234","SOL235","SOL236","SOL237","SOL238","SOL239","SOL240","SOL241","SOL242","SOL243","SOL244","SOL245","SOL246","SOL247","SOL248","SOL249","SOL250","SOL251","SOL252","SOL253","SOL254","SOL255","SOL256","SOL257","SOL258","SOL259","SOL260","SOL261","SOL262","SOL263","SOL264","SOL265","SOL266","SOL267","SOL268","SOL269","SOL270","SOL271","SOL272","SOL273","SOL274","SOL275","SOL276","SOL277","SOL278","SOL279","SOL280","SOL281","SOL282","SOL283","SOL284","SOL285","SOL286","SOL287","SOL288","SOL289","SOL290","SOL291","SOL292","SOL293","SOL294","SOL295","SOL296","SOL297","SOL298","SOL299","SOL300","SOL301","SOL302","SOL303","SOL304","SOL305","SOL306","SOL307","SOL308","SOL309","SOL310","SOL311","SOL312","SOL313","SOL314","SOL315","SOL316","SOL317","SOL318","SOL319","SOL320","SOL321","SOL322","SOL323","SOL324","SOL325","SOL326","SOL327","SOL328","SOL329","SOL330","SOL331","SOL332","SOL333","SOL334","SOL335","SOL336","SOL337","SOL338","SOL339","SOL340","SOL341","SOL342","SOL343","SOL344","SOL345","SOL346","SOL347","SOL348","SOL349","SOL350","SOL351","SOL352","SOL353","SOL354","SOL355","SOL356","SOL357","SOL358","SOL359","SOL360","SOL361","SOL362","SOL363","SOL364","SOL365","SOL366","SOL367","SOL368","SOL369","SOL370","SOL371"];
    
    %371 molecules composition below (evenly concentrated case)
    sysInfo.yf = xlsread('dataConcentration');
    
    %Albatross-5-12comp feed composition below
    %sysInfo.yf = [0.9;0.1];  % composition spec (must sum to 1)
    
    %Albatross-5-12comp feed composition below
    %sysInfo.yf = [0.0193;0.8924;0.0064;0.0092;0.0077;0.0117;0.0107;0.0071;0.0084;0.0074;0.0098;0.01];  % composition spec (must sum to 1)
    
    %Ducky-10-12comp feed composition below
    %sysInfo.yf = [0.0126767015392605;0.874955024936112;0.0107799133222996;0.0112794710793953;0.0109978527000903;0.0111994143059899;0.011689929770511;0.0111221962999486;0.0111741546080492;0.0109983150008092;0.0114540441780076;0.0116729822595269];  % composition spec (must sum to 1)
    
    %Ducky-9-12comp feed composition below
    %sysInfo.yf = [0.0103276162247238;0.883434986215907;0.0102433116826168;0.0104430839738515;0.0101305616107897;0.0109774377992201;0.011252857686684;0.0107043883954339;0.0104608400849019;0.00927498940308899;0.011060428809073;0.0116894981137093];  % composition spec (must sum to 1)
    
    %Mallard-1-12comp feed composition below
    %sysInfo.yf = [0.01370891122;0.872133954; 0.01169827194;0.01152161375;0.0115658596;0.01264385258;0.01316118212;0.00509529692;0.01210352992;0.01168821572;0.01196063483;0.0127186774];  % composition spec (must sum to 1)
    
    PuM = 30; %feed-side pressure (bar)
  %system specs
    sysInfo.n = length(sysInfo.yf);  % number of permeants
    sysInfo.Pu = PuM*0.9869;  % feed side pressure [bar]*0.9832 = [atm]
    sysInfo.Pd = 1;  % support side pressure (atm)
    sysInfo.pervapMode = 0; % %BETA(UNSTABLE) pervap capability -- yes if == 1 (assume downstream pressure = 0 bar)
    sysInfo.T = 295;  % system temperature (K)
    sysInfo.R = 82.05;  % gas constant (atm cm^3/ mol K)
    sysInfo.diffFit = 0; %see code below for diffusivity fitting
    sysInfo.crossDiffFudge = 0; %if = 0 then no specified cross-diffusivities D_ij will be used
    if sysInfo.crossDiffFudge == 1
        sysInfo.crossDiffSpecs = [1,2,1,3]; % e.g. [i,j] = D_ij^MS
        sysInfo.crossDiffVals = [0.1,0.5]; %um2/s
    end
    sysInfo.lowDiffErrorBar = 0; % if = 1 then diffusivities fit to low error sing comp will be used
%------------------------------------------------------------------------------------------------------------------------------------%
end
