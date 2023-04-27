%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemDiffFitSolve(sysInfo,diff,singFlux)                  %
% Description: Solve for singCompDiffs given singCompFluxes               %
% Input:       sysInfo            - struct defining simulation specs      %         
%                                (see asyMemLocalSolve function for specs)%
%              singFluxVec        - experimental flux values to fit       %
%                                     singCompDiff [LMH]                  %
%              singFluxVec_Error  - experimental flux error values to fit %
%                                     singCompDiff [LMH]                  %
%              modSpec            - vector of number values to choose fit %
%                                     sorp model (1-FH, 2-DSM, 3-FH-LM)   %
% Output:      allModSingCompDiff - corresponding diff values [um^2/s]    % 
%              allModSingCompDiff - corresponding LOW diff values[um^2/s] % 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [allModSingCompDiff,allModSingCompDiffLOW] = asyMemDiffFitSolve(sysInfo,singFluxVec,singFluxVec_Error,modSpec)
    sysInfo.n = 1;
    sysInfo.yf = [1];
    singFluxVec_LOW = [singFluxVec-singFluxVec_Error]; %LMH
    sysInfo.diffFit = 1;
    allModSingCompDiff = [];
    allModSingCompDiffLOW = [];
    for i = modSpec %change for how many models want fit (i.e i = 1:2 for FH, DSM, and FH-LM)
        if i == 1
            sysInfo.memPhaseModel = 'F-H';
        elseif i == 2
            sysInfo.memPhaseModel = 'DSM';
        else
            sysInfo.memPhaseModel = 'FH-LM';
        end
        singCompDiff = [];
        singCompDiffLOW = [];
        for j = 1:length(singFluxVec)
            sysInfo.mixID = sysInfo.mixID_OG(j);
            [singCompDiff_j] = fzero(@(diff)asyMemDiffFit_RHS(sysInfo,diff,singFluxVec(j)),1);
            [singCompDiffLOW_j] = fzero(@(diff)asyMemDiffFit_RHS(sysInfo,diff,singFluxVec_LOW(j)),1);
            singCompDiff = [singCompDiff;singCompDiff_j];
            singCompDiffLOW = [singCompDiffLOW;singCompDiffLOW_j];
        end
        allModSingCompDiff = [allModSingCompDiff;singCompDiff];
        allModSingCompDiffLOW = [allModSingCompDiffLOW;singCompDiffLOW];
    end
end

