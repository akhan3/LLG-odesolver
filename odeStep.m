function [Mnext] = odeStep(sp, bc, M, Hext)
    if sp.odeSolver_rk4
        Mnext = rk4Step(sp, bc, M, Hext);
    else 
        Mnext = eulerStep(sp, bc, M, Hext);
    end
end


function [Mnext] = eulerStep(sp, bc, M, Hext)
    H = Hfield(sp, bc, M, Hext);
    Mprime = LLG_Mprime(sp, M, H);
    Mnext = M + Mprime*sp.dt;
end


function [Mnext] = rk4Step(sp, bc, M, Hext)
    % k1
    H = Hfield(sp, bc, M, Hext);
    k1 = LLG_Mprime(sp, M, H);
    Mnext = M + k1*sp.dt/2.0;
    % k2
    H = Hfield(sp, bc, Mnext, Hext);    % TODO: error bcoz not interpolating Hext(t+dt/2)
    k2 = LLG_Mprime(sp, Mnext, H);
    Mnext = M + k2*sp.dt/2.0;
    % k3
    H = Hfield(sp, bc, Mnext, Hext);    % TODO: error bcoz not interpolating Hext(t+dt/2)
    k3 = LLG_Mprime(sp, Mnext, H);
    Mnext = M + k3*sp.dt;
    % k4
    H = Hfield(sp, bc, Mnext, Hext);    % TODO: error bcoz not using Hext(t+dt)
    k4 = LLG_Mprime(sp, Mnext, H);
    Mnext = M + (k1+2*k2+2*k3+k4)/6.0 * sp.dt;
end


function [H] = Hfield(sp, bc, M, Hext)
    H = zeros(size(M));   % pre-allocate Hfield to tell the size
    if sp.useGPU
        H = gpuArray(H);
    end
    %% Coupling field for boundary elements
    H(:,1  ,:) = sp.c * bc.Mtop(:,:);   % top row
    H(:,end,:) = sp.c * bc.Mbot(:,:);   % bottom row
    H(:,:,end) = sp.c * bc.Mrig(:,:);   % right column
    H(:,:  ,1) = sp.c * bc.Mlef(:,:);   % left column

    %% Coupling field for rest of the matrix
    Mtop = M(:,1:end-2,2:end-1);
    Mbot = M(:,3:end,  2:end-1);
    Mrig = M(:,2:end-1,1:end-2);
    Mlef = M(:,2:end-1,3:end  );
    H(:,2:end-1,2:end-1) = sp.c * (Mtop + Mbot + Mrig + Mlef);

    %% Interaction with itself (Demagnetizating field)
    H(1,:,:) = H(1,:,:) + -sp.N(1,1) * M(1,:,:);
    H(2,:,:) = H(2,:,:) + -sp.N(2,2) * M(2,:,:);
    H(3,:,:) = H(3,:,:) + -sp.N(3,3) * M(3,:,:);
end


function [Mprime] = LLG_Mprime(sp, M, H)
%     McrossH = cross(M,H);
%     Mprime = -sp.gamma * McrossH - (sp.alpha*sp.gamma/sp.Ms) * cross(M,McrossH);
    [MH(1,:,:), MH(2,:,:), MH(3,:,:)]= crossProduct(M,H);
    [MMH(1,:,:), MMH(2,:,:), MMH(3,:,:)]= crossProduct(M,MH);
    Mprime = -sp.gamma * MH - (sp.alpha*sp.gamma/sp.Ms) * MMH;
end


function [C1, C2, C3] = crossProduct(A,B)
    C1 = A(2,:,:).*B(3,:,:) - A(3,:,:).*B(2,:,:);
    C2 = A(3,:,:).*B(1,:,:) - A(1,:,:).*B(3,:,:);
    C3 = A(1,:,:).*B(2,:,:) - A(2,:,:).*B(1,:,:);
end
