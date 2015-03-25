function [pcCoeff, spatialBase] = makePcMeanFree(pcCoeffIn, spatialBaseIn)
%Usually, the first stiffness matrix calculated from the NORTA expansion of a pde coefficient field corresponds to the mean (the first KL eigenfunction). This is not the case if one does not use the KL expansion, for example when expanding a piecewise constant random function into random variables and spatial indicator functions. 
%The topic is of importance because Elmar zanders preconditioning substracts the stiffness matrix corresponding to the mean from the other matrices and assumes it to be the first block in the stochastic operator.
%This function serves that case by constructing a "mean variable" variable separated into a spatial vector of means and a degree-zero PC expansion, and putting the first in first column of spatial functions, and for the second add a column and a line to the multiindex matrix and a one into the upper left corner.
 nH_alpha = size(pcCoeffIn,2)
 nRandomVars = size(pcCoeffIn,1)+1

means = pcCoeffIn(:,1);

 pcCoeff = (zeros(nRandomVars,nH_alpha))
 pcCoeff(1,1) =1; 
 pcCoeff(2:nRandomVars,2:nH_alpha) =pcCoeffIn(:,2:nH_alpha); 
%  for i_randVar = 2:nRandomVars
% 	pc
% end
 
 
 nPts =size(spatialBaseIn,1)
 spatialBase = zeros(nPts, nRandomVars);
%debug: pcBase(:, 1)= pcBaseIn(:,1)~=0;

%spatialBase(:, 1)= means(1)*(spatialBaseIn(:,1)~=0);
for i=1:nRandomVars-1
    spatialBase(:, 1)= spatialBase(:, 1) + means(i)*(spatialBaseIn(:,i)~=0);
end
%i.e. domainwise means are written to first spatial basis vector
spatialBase(:, 2:nRandomVars)=spatialBaseIn;