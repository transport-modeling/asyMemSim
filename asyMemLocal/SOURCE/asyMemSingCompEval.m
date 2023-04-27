%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemSingCompEval(sysInfoExt,pressureVec,modSpec)         %
% Description: Evaluate singComp mixture for modSpecs and PressureVec     %
% Input:       sysInfo     - struct defining simulation specs             %      
%                              (see asyMemLocalSolve function for specs)  %
%              pressureVec - vecotr of pressures to run through           %
%              modSpec     - vector of number values to choose fit        %
%                              sorp model (1-FH, 2-DSM, 3-FH-LM)          %
% Output:      allMPSCP    - matrix of singComp fluxes for each pressure  %
%                              and each model                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allMPSCF,failedRuns] = asyMemSingCompEval(sysInfo,pressureVec,modSpec)
    allMPSCF = [];
    failedRuns = [];
    for k = 1:length(pressureVec)
        sysInfo.Pu = pressureVec(k)*0.9869;
        allModSingCompFlux = [];
        for i = modSpec %change for how many models want fit (i.e i = 1:2 for FH, DSM, and FH-LM)
            if i == 1
                sysInfo.memPhaseModel = 'F-H';
            elseif i == 2
                sysInfo.memPhaseModel = 'DSM';
            else
                sysInfo.memPhaseModel = 'FH-LM';
            end
            singCompFlux = [];
            for j = 1:length(sysInfo.mixID_OG)
                sysInfo.mixID = sysInfo.mixID_OG(j);
                [localCompFlux,partialFlux,exitFlag] = asyMemLocalSolve(sysInfo);
                singCompFlux = [singCompFlux;partialFlux];
                if exitFlag ~= 1 && exitFlag ~= 2 && exitFlag ~= 3 && exitFlag ~= 4
                    failedRuns = [failedRuns,[exitFlag;k;i]];
                end
            end
            allModSingCompFlux = [allModSingCompFlux;singCompFlux];
        end
       allMPSCF = [allMPSCF,allModSingCompFlux];
    end
end



