function F = steady_ntwsoe_system(x,alpha_vec, alphaV_vec, beta_mat, MCi, PMi, PL, PV, CHi, Xi,epsY_vec, epsM_vec, A_vec,pHvec)
    % Unpack variables from x
    % Example: x = [M; L; Vi; Yi; ...] for all sectors
    nsec = length(alpha_vec);
    M   = x(1:nsec);
    L   = x(nsec+1:2*nsec);
    Vi  = x(2*nsec+1:3*nsec);
    Yi  = x(3*nsec+1:4*nsec);
    
    % Compute intermediate use correctly: sum_j beta(j,i) * (PM(j)/PH(i))^epsM(j) * M(j)
    intermediate_use = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            intermediate_use(i) = intermediate_use(i) + beta_mat(j,i) * (PMi(j)/pHvec(i))^epsM_vec(j) * M(j);
        end
    end

    F1 = M  -  (MCi./PMi).^epsY_vec.*alpha_vec.*(CHi + Xi + intermediate_use);

    F2 = L - (MCi./PL).^epsY_vec.*(1-alpha_vec-alphaV_vec).* ...
        (CHi + Xi + intermediate_use);
    F3 = Vi - (MCi./PV).^epsY_vec.*alphaV_vec.* ...
        (CHi + Xi + intermediate_use);
    F4 = Yi - A_vec.*(alpha_vec.^(1./epsY_vec).*M.^((epsY_vec-1)./epsY_vec) ...
        + alphaV_vec.^(1./epsY_vec).*Vi.^((epsY_vec-1)./epsY_vec) ...
        + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L.^((epsY_vec-1)./epsY_vec)).^(epsY_vec./(epsY_vec-1));
    % Add more equations as needed for your system

    % Concatenate all equations
    F = [F1; F2; F3; F4];
end