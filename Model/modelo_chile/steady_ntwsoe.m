
function residual=steady_ntwsoe(x_vec,PVstar,epsilon,...
    varrho_vec,sigmaH,gammag_vec,gammas_vec,om_g,om_s,chiX_vec,omegaX,etastar,Ystar,...
    alpha_vec,alphaV_vec,beta_mat,epsY_vec,epsM_vec,GAMMA,CHI,PSI,A_vec,tb_target)
    nsec=length(varrho_vec);
    pHvec = x_vec(1:nsec,1);
    w     = x_vec(nsec+1,1);
    q     = x_vec(nsec+2,1);
    C     = x_vec(nsec+3,1);

   
    PL = ones(nsec,1)*w; %fine
    PV=q*PVstar; %fine
    
    MCi=(epsilon-1)/epsilon*pHvec; % Marginal costs given prices  %fine
    PMi = (beta_mat(:,:)*pHvec.^(1-epsM_vec)).^(1./(1-epsM_vec)); % Price index of intermediates by sector.  %fine
    
    pvec=(varrho_vec.^(sigmaH).*pHvec.^(1-sigmaH)+(1-varrho_vec).^(sigmaH).*PV.^(1-sigmaH)).^(1/(1-sigmaH)); % price of the good by the side of consumers

    p_g = prod(pvec.^gammag_vec); %multiply Goods price index
    p_s = prod(pvec.^gammas_vec); %multiply Services price index

    C_g = om_g*C/p_g; %consumption of goods
    C_s = om_s*C/p_s; % consumption of services

    C_gi = gammag_vec.*(p_g./pvec).*C_g; % consumption of each good i if it is in goods' set
    C_si = gammas_vec.*(p_s./pvec).*C_s; % consumption of each good i if it is in services' set

    CH_gi = varrho_vec.^(sigmaH).*(pHvec./pvec).^(-sigmaH).*C_gi; % consumption of good i in goods if home produced
    CH_si = varrho_vec.^(sigmaH).*(pHvec./pvec).^(-sigmaH).*C_si; % consumption of good i in services if home produced

    CF_gi = (1-varrho_vec).^(sigmaH).*(PV./pvec).^(-sigmaH).*C_gi;  % consumption of good i in goods if foreign produced
    CF_si = (1-varrho_vec).^(sigmaH).*(PV./pvec).^(-sigmaH).*C_si;  % consumption of good i in services if foreign produced

    CHi = CH_gi+CH_si; % home demand for good i
    CFi = CF_gi+CF_si; % foreign demand for good i

    PX =prod(pHvec.^chiX_vec); % Exports price index

    X = omegaX*(PX/q)^(-etastar)*Ystar; % Total exports

    Xi  = chiX_vec.*X.*PX./pHvec; % sector i exports

   % Compute intermediate use correctly for initial guess
    intermediate_use_guess = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            intermediate_use_guess(i) = intermediate_use_guess(i) + beta_mat(j,i);
        end
    end
    intermediate_use_guess = intermediate_use_guess .* mean(CHi + Xi);

    M   = (MCi./PMi).^epsY_vec.*alpha_vec.*(CHi + Xi + intermediate_use_guess);
    
    % Compute proper intermediate use with M
    intermediate_use = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            intermediate_use(i) = intermediate_use(i) + beta_mat(j,i) * (PMi(j)/pHvec(i))^epsM_vec(j) * M(j);
        end
    end
    
    L   = (MCi./PL).^epsY_vec.*(1-alpha_vec-alphaV_vec).*(CHi + Xi + intermediate_use);
    Vi   = (MCi./PV).^epsY_vec.*alphaV_vec.*(CHi + Xi + intermediate_use);
    Yi  = A_vec.*(alpha_vec.^(1./epsY_vec).*M.^((epsY_vec-1)./epsY_vec)... 
        + alphaV_vec.^(1./epsY_vec).*Vi.^((epsY_vec-1)./epsY_vec)...
        + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L.^((epsY_vec-1)./epsY_vec)).^(epsY_vec./(epsY_vec-1));
    x_guess1=[M; L; Vi; Yi];
    options1 = optimoptions('fsolve','TolFun',1e-10,'Display','off');

    x_sol=fsolve(@(x_vec) steady_ntwsoe_system(x_vec,alpha_vec, alphaV_vec, beta_mat, MCi, PMi, PL, PV, CHi, Xi,epsY_vec, epsM_vec, A_vec,pHvec),x_guess1,options1);
    M   = x_sol(1:nsec);
    L   = x_sol(nsec+1:2*nsec);
    Vi  = x_sol(2*nsec+1:3*nsec);
    Yi  = x_sol(3*nsec+1:4*nsec);


    %TCR = P*/CPI
    TB = PX*X-PV*(sum(CFi)+sum(Vi)); %PX*x_gdp-PV*V/GDP
    GDP=C+TB;
    N=(C^(-GAMMA)*w/CHI)^(1/PSI);

    % Compute intermediate use for market clearing residual
    intermediate_use_final = zeros(nsec,1);
    for i = 1:nsec
        for j = 1:nsec
            intermediate_use_final(i) = intermediate_use_final(i) + beta_mat(j,i) * (PMi(j)/pHvec(i))^epsM_vec(j) * M(j);
        end
    end

    residual=zeros(nsec+3,1);
    residual(1:nsec,1) = Yi-CHi - Xi - intermediate_use_final; %this is nx1 vector Prices solve this    
    residual(nsec+1,1) = TB/GDP-tb_target; %consumption?
    residual(nsec+2,1) = N-sum(L); %wages
    %residual(nsec+3,1) = 1-p_g^om_g*p_s^om_s; %real exchange rate
    residual(nsec+3,1) = C -(C_g/om_g)^om_g*(C_s/om_s)^om_s; %consumption definition    