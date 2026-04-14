function resid=pricesmc(pvec,MC,PL,epsilon,beta_mat,epsM_vec,A_vec,alpha_vec,epsY_vec)

    for j=1:length(pvec)
        PM(j) = sum(beta_mat(j,:)*pvec.^(1-epsM_vec(j)))^(1/(1-epsM_vec(j))); 
        MC(j) = 1/A_vec(j)*(alpha_vec(j)*PM(j)^(1-epsY_vec(j))+(1-alpha_vec(j))*PL(j)^(1-epsY_vec(j)))^(1/(1-epsY_vec(j)));    
    end

    resid = pvec - epsilon/(epsilon-1)*MC;
