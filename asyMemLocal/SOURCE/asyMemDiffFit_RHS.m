%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function:    asyMemDiffFit_RHS(sysInfoExt,diff,singFlux)                %
% Description: Evaluate RHS function for solving singCompDiffs            %
% Input:       sysInfoExt - external struct defining simulation specs     %
%                             (see asyMemLocalSolve function for specs)   %
%              diff       - iteration diffusivity [um^2/s]                %
%              singFlux   - experimental flux values to fit singCompDiff  %
% Output:      funRHS     - function values vector for diff fit equations %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [funRHS] = asyMemDiffFit_RHS(sysInfo,diff,singFlux)
    sysInfo.diffFitIter = diff;
    localCompFlux = asyMemLocalSolve(sysInfo);
    funRHS = singFlux - localCompFlux(2);
end

