function [Mnext] = odeStep(sp, bc, M, Hext)

% [Mnext] = odeStepMatlab(sp, bc, M, Hext);

%MM = single(reshape(1:numel(M), size(M)) + .1);
%HHext = single(reshape(1:numel(M), size(M)) + 100.1);
%MM(1,:) = 100+MM(1,:);
%MM(2,:) = 200+MM(2,:);
%MM(3,:) = 300+MM(3,:);
%squeeze(MM(1,:,:));
%MM

[Mnext] = odeStepComp(sp, single(M), single(Hext));
%[Mnext] = odeStepMat(sp, bc, single(M), single(Hext));



end
