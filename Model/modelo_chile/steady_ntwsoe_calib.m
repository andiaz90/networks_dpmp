function residual=steady_ntwsoe_calib(x_vec,PVstar,epsilon,...
    sigmaH,gammag_vec,gammas_vec,om_g,om_s,chiX_vec,etastar,Ystar,...
    alpha_vec,alphaV_vec,beta_mat,epsY_vec,epsM_vec,GAMMA,CHI,PSI,...
    A_vec,tb_target,CHshare_target,X_share_target)
    
    nsec=length(chiX_vec);

    pHvec      = x_vec(1:nsec,1);
    w          = x_vec(nsec+1,1);
    q          = x_vec(nsec+2,1);
    C          = x_vec(nsec+3,1);
    varrho_vec = x_vec(nsec+4:2*nsec+3,1);
    omegaX     = x_vec(2*nsec+4,1);

   
    PL = ones(nsec,1)*w;
    PV=q*PVstar;
    
    MCi=(epsilon-1)/epsilon*pHvec; % Marginal costs given prices
    PMi = (beta_mat(:,:)*pHvec.^(1-epsM_vec)).^(1./(1-epsM_vec)); % Price index of intermediates by sector. 
    
    pvec=(varrho_vec.*pHvec.^(1-sigmaH)+(1-varrho_vec).*PV.^(1-sigmaH)).^(1/(1-sigmaH)); % price of the good by the side of consumers

    p_g = prod(pvec.^gammag_vec); %multiply Goods price index
    p_s = prod(pvec.^gammas_vec); %multiply Services price index

    C_g = om_g*C/p_g; %consumption of goods
    C_s = om_s*C/p_s; % consumption of services

    C_gi = gammag_vec.*(p_g./pvec).*C_g; % consumption of each good i if it is in goods' set
    C_si = gammas_vec.*(p_s./pvec).*C_s; % consumption of each good i if it is in services' set

    CH_gi = varrho_vec.*(pvec./pHvec).*C_gi; % consumption of good i in goods if home produced
    CH_si = varrho_vec.*(pvec./pHvec).*C_si; % consumption of good i in services if home produced

    CF_gi = (1-varrho_vec).*(pvec./PV).*C_gi;  % consumption of good i in goods if foreign produced
    CF_si = (1-varrho_vec).*(pvec./PV).*C_si;  % consumption of good i in services if foreign produced

    CHi = CH_gi+CH_si; % home demand for good i
    CFi = CF_gi+CF_si; % foreign demand for good i

    PX =prod(pHvec.^chiX_vec); % Exports price index

    X = omegaX*(PX/q)^(-etastar)*Ystar; % Total exports

    Xi  = chiX_vec.*X.*PX./pvec; % sector i exports

    mkupV_ss = epsilonV./(epsilonV-1);
   % next, we need to obtain the inputs demanded by each sector. We assume
   % the demand is "right" and use it as the levet of production to
   % generate. Sthill, this is a nonlinear problem because Y depends on M
   % and M depends on Y nonlinearly
    M   = 1./(1-alpha_vec.*beta_mat(:,:)*(MCi./PMi).^epsY_vec.*(PMi./pHvec).^epsM_vec).*(MCi./PMi).^epsY_vec.*alpha_vec.*(CHi + Xi);
    DDi = CHi + Xi +beta_mat(:,:)*M.*(PMi./pHvec).^epsM_vec;
    L   = (MCi./PL).^epsY_vec.*(1-alpha_vec-alphaV_vec).*DDi;
    Vi   = (MCi./PV).^epsY_vec.*alphaV_vec.*DDi;


    Yi  = A_vec.*(alpha_vec.^(1./epsY_vec).*M.^((epsY_vec-1)./epsY_vec)... 
        + alphaV_vec.^(1./epsY_vec).*Vi.^((epsY_vec-1)./epsY_vec)...
        + (1-alphaV_vec-alpha_vec).^(1./epsY_vec).*L.^((epsY_vec-1)./epsY_vec)).^(epsY_vec./(epsY_vec-1));

    %TCR = P*/CPI
    TB = PX*X-PV*(sum(CFi)+sum(Vi)); %PX*x_gdp-PV*V/GDP
    GDP=C+TB;
    N=(C^(-GAMMA)*w/CHI)^(1/PSI);

    CHshare = pHvec.*CHi./(pvec.*C_gi+pvec.*C_si);

    X_share = PX*X/GDP;

    residual=zeros(nsec+3+nsec,1);
    residual(1:nsec,1) = Yi-DDi; %this is nx1 vector Prices solve this    
    residual(nsec+1,1) = TB/GDP-tb_target; %consumption?
    residual(nsec+2,1) = N-sum(L); %wages
    residual(nsec+3,1) = 1-p_g^om_g*p_s^om_s; %real exchange rate
    residual(nsec+4:end,1) = CHshare-CHshare_target; %real exchange rate
    
    residual(nsec+4:end,1) = X_share-X_share_target; %real exchange rate
