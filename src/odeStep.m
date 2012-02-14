function [Mnext] = odeStep(sp, bc, M, Hext)

    [Mnext] = odeStepComp(sp, single(M), single(Hext));

end
