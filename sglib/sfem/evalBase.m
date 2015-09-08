%evaluates arbitrary spatial trigonometric base.
%Caller has to take care that X and Y are scaled and shifted to interval
%[-1,1] -this function has no chance to do that for arbitrary points!
%TODO: Think of giving interval of period as arg!
function [spatialBase]=evalBase(k, phase, X, Y)
for ik=1:size(k,1)
    for jk=1:size(k,2) %reshape(k(ik,jk,1:2),1,[])*[X;Y]
     spatialBase(ik,jk,:)=cos( pi*( k(ik,jk,1)*X+phase(ik,jk,1) + k(ik,jk,2)*Y+ phase(ik,jk,2)));%[X';Y']+phase(ik,jk,1:2)));
    end
end
end
