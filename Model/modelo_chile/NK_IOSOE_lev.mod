//% N sector NK model with Input-Output Linkages

// Include blocks defined elsewhere
@#include "definition_block_nsec.mod" 
@#include "definition_block_io.mod" 
@#include "definition_block_lab.mod" 

// Define variables
var chi Lab_costs Price_costs w VA C Ctot Ctotg Ctots M C_g C_s p_g p_s N pi pi_g pi_s r om_g om_s vi Y
PV r_star Rworld pi_e Bstar Pipstar Q TB PX X Ystar xi_psi V PVstar GDP IMP mkupV CFs CFg

    @#for i in 1:nsec
        PH_@{i}
    @#endfor
    @#for i in 1:nsec
        P_@{i}
    @#endfor
    
    @#for i in 1:nsec
        L_@{i}
    @#endfor
    @#for i in 1:nsec
        M_@{i}
    @#endfor
    @#for i in 1:nsec
        Y_@{i}
    @#endfor
    @#for i in 1:nsec
        PM_@{i}
    @#endfor
    @#for i in 1:nsec
        PL_@{i}
    @#endfor
    @#for i in 1:nsec
        MC_@{i}
    @#endfor
    @#for i in 1:nsec
        A_@{i}
    @#endfor
    @#for i in 1:nsec
        V_@{i}
    @#endfor
    @#for i in 1:nsec
    CHg_@{i} 
    CFg_@{i} 
    CHs_@{i} 
    CFs_@{i} 
    Cg_@{i} 
    Cs_@{i}
    @#endfor
    ;

varexo eps_om eps_i epschi eps_psi
    @#for i in 1:nsec
        epsA_@{i}
    @#endfor
;

parameters gamma psi beta phi epsilon rho rho_om1 rho_tfp1 rho_om2 rho_tfp2 ombar
rhoi rhoirule ilabcosts sigma_i sigma_om sigma_L_agg Rworld_ss
bbar chii_b omegaX epsilonX ystar_ss etastar kappaV epsilonV sigmaH Pistar_ss PVstar_ss
IMP_ss rho_psi sigma_psi

    @#for i in 1:nsec
        gammag_@{i}
        alpha_@{i}
        gammas_@{i}
        @#for j in 1:nsec
            beta_@{i}_@{j}
        @#endfor
        epsY_@{i}
        epsM_@{i}
        kappa_@{i}
        dummyg_@{i}
        dummys_@{i}
        cl_@{i}
        clneg_@{i}
        cm_@{i}
        chiX_@{i}
        alphaV_@{i}
        varrho_@{i}
        isigma_tfp_@{i}
        PL_ss@{i}
    @#endfor
    ;  

      
load params_val.mat;
//load params_val_SOE.mat;

sigma_i     = sigma_i_val;
sigma_L_agg = sigma_L_agg_val;
sigma_om    = sigma_om_val;
ilabcosts   = ilabcosts_val;

gamma      = gamma_val;
psi        = 1;
beta       = beta_val;
phi        = phi_val;
epsilon    = 10;
rho        = rho_val;
rho_om1    = rho_om1_val;
rho_om2    = rho_om2_val;
rho_tfp1   = rho_tfp1_val;
rho_tfp2   = rho_tfp2_val;
rhoi       = 0.5+0*rhoi_val;
rhoirule   = rhoirule_val;
ombar      = ombar_val;
chii_b     = chii_b_val;
bbar       = bbar_val;
epsilonX   = epsilonX_val;
omegaX     = omegaX_val;
ystar_ss   = ystar_ss_val;
etastar    = etastar_val;
epsilonV   = epsilonV_val;
kappaV     = kappaV_val;
sigmaH     = sigmaH_val;
w_ss       = w_ss;
C_ss       = C_ss;
GDP_ss     = GDP_ss;
N_ss       = N_ss;
p_s_ss     = p_s_ss;
p_g_ss     = p_g_ss;
C_s_ss     = C_s_ss;
C_g_ss     = C_g_ss;
Pistar_ss  = Pistar_ss_val;
PVstar_ss  = PVstar_ss_val;
Rworld_ss  = Rworld_ss_val;
Bstar_ss   = Bstar_ss;
IMP_ss     = IMP_ss_val;
rho_psi    = rho_psi_val;
sigma_psi  = sigma_psi_val;
@#for i in 1:nsec
    gammag_@{i} = modgammag(@{i});
    gammas_@{i} = modgammas(@{i});
    alpha_@{i}  = modalpha(@{i});
    @#for j in 1:nsec
        beta_@{i}_@{j} = modbeta(@{i},@{j});
    @#endfor
    epsY_@{i}   = modepsY(@{i});
    epsM_@{i}   = modepsM(@{i});
    kappa_@{i}  = modkappa(@{i});
    dummyg_@{i} = goods(@{i});
    dummys_@{i} = services(@{i});
    cl_@{i}     = modcl(@{i});
    clneg_@{i}  = modclneg(@{i});
    cm_@{i}     = modcm(@{i});
    chiX_@{i}   = modchiX(@{i});
    alphaV_@{i} = modalphaV(@{i});
    varrho_@{i} = modvarrho(@{i});
    isigma_tfp_@{i} = isigma_tfp_val(@{i});
    PL_ss@{i}   = PL_ss(@{i});
    CFg_ss@{i}  = CFg_ss(@{i});
    CFs_ss@{i}  = CFs_ss(@{i});
    CHg_ss@{i}  = CHg_ss(@{i});
    CHs_ss@{i}  = CHs_ss(@{i});
    Vi_ss@{i}   = Vi_ss(@{i});
    PH_ss@{i}   = pH_ss(@{i});
    MC_ss@{i}   = MCi_ss(@{i});
    Y_ss@{i}   = Yi_ss(@{i});
    L_ss@{i}   = L_ss(@{i});
    Cgi_ss@{i} = C_gi_ss(@{i});
    Csi_ss@{i} = C_si_ss(@{i});
    P_ss@{i}   = P_ss(@{i});
    PMi_ss@{i} = PMi_ss(@{i});
    Mi_ss@{i}  = M_ss(@{i});
@#endfor

model(bytecode);

//% Overall Consumption Basket

//% Goods/services aggregator
C = (C_g/exp(om_g))^exp(om_g)*(C_s/exp(om_s))^exp(om_s) ;

//% Goods/services demands
p_g = exp(om_g)*C/C_g;
p_s = exp(om_s)*C/C_s;


% Good j price index.

@#for j in 1:nsec
    P_@{j}=(varrho_@{j}^(sigmaH)*PH_@{j}^(1-sigmaH)+(1-varrho_@{j})^(sigmaH)*PV^(1-sigmaH))^(1/(1-sigmaH));
@#endfor

//% Relative Price of Goods
p_g = (1
@#for j in 1:nsec
    *(P_@{j})^gammag_@{j}
@#endfor
);

//% Relative Price of Services
p_s = (1
@#for j in 1:nsec
    *(P_@{j})^gammas_@{j}
@#endfor
);
@#for j in 1:nsec
    Cg_@{j}=gammag_@{j}*(p_g/P_@{j})*C_g;
    Cs_@{j}=gammas_@{j}*(p_s/P_@{j})*C_s;
@#endfor
% Demands for goods
@#for j in 1:nsec   
    CHg_@{j} =    varrho_@{j}^(sigmaH)*(PH_@{j}/P_@{j})^(-sigmaH)*Cg_@{j};
    CFg_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV/P_@{j})^(-sigmaH)*Cg_@{j};    
    CHs_@{j} =    varrho_@{j}^(sigmaH)*(PH_@{j}/P_@{j})^(-sigmaH)*Cs_@{j};
    CFs_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV/P_@{j})^(-sigmaH)*Cs_@{j};
@#endfor

//% Household Euler Equation
C^(-gamma) = beta*(C(+1)^-gamma)*r/pi(+1);

//% Relative Inflation: Goods vs Services
pi_g*C_g/C_g(-1)*(exp(om_s)/exp(om_s(-1))) = pi_s*C_s/C_s(-1)*(exp(om_g)/exp(om_g(-1)));

//% Labor Market Clearing
N = (
    @#for i in 1:nsec
        +L_@{i} 
        +cl_@{i}/2*L_@{i}*(L_@{i}/L_@{i}(-1)-1)^2
    @#endfor
    );

Lab_costs = (
    @#for i in 1:nsec
        + ilabcosts*cl_@{i}/2*L_@{i}*(L_@{i}/L_@{i}(-1)-1)^2
    @#endfor
    );
    
//% Labor leisure condition
C^(-gamma)*w = chi*N^psi;


//% Relative Inflation: Goods vs Overall
pi_g*C_g/C_g(-1) = exp(om_g)/exp(om_g(-1))*pi*C/C(-1);

//% Gross Output
Y = (
    @#for i in 1:nsec
        +Y_@{i}
    @#endfor
    );

//% Total Consumption (Note: = Value Added)
Ctot = (
@#for i in 1:nsec
    +gammag_@{i}*C_g*p_g/P_@{i} + gammas_@{i}*C_s*p_s/P_@{i}
@#endfor
);

//% Total Consumption: Services
Ctots = (
@#for i in 1:nsec
    +gammas_@{i}*C_s*p_s/P_@{i}
@#endfor
);

//% Total Consumption: Goods
Ctotg = (
@#for i in 1:nsec
    +gammag_@{i}*C_g*p_g/P_@{i}
@#endfor
);

    //% Value Added
    VA = (
        @#for i in 1:nsec
            +Y_@{i} - M_@{i}
        @#endfor
    );

GDP = C + TB;
//% Total Intermediate Use
M = (
@#for i in 1:nsec
    @#for j in 1:nsec
        + beta_@{j}_@{i}*(PM_@{j}/PH_@{i})^epsM_@{j}*M_@{j}
    @#endfor
@#endfor
);

//% Total Price Adjustment Costs
Price_costs = (
@#for i in 1:nsec
        + (kappa_@{i}/2)*(pi*PH_@{i}/PH_@{i}(-1)-1)^2*Y_@{i}
@#endfor
);


//% Sector-by-Sector Equations:
@#for i in 1:nsec

    //% Intermediates Price Index
    PM_@{i} = (
    @#for j in 1:nsec
        + beta_@{i}_@{j}*(PH_@{j})^(1-epsM_@{i})
    @#endfor
    )^(1/(1-epsM_@{i})); 

    //% Production Function
        Y_@{i} = exp(A_@{i})*((alpha_@{i})^(1/epsY_@{i})*(M_@{i})^((epsY_@{i}-1)/epsY_@{i}) 
                               + (alphaV_@{i})^(1/epsY_@{i})*(V_@{i})^((epsY_@{i}-1)/epsY_@{i})
                               + (1-alpha_@{i}-alphaV_@{i})^(1/epsY_@{i})*(L_@{i})^((epsY_@{i}-1)/epsY_@{i}))^(epsY_@{i}/(epsY_@{i}-1));

    //% Intermediates Demand
        MC_@{i}*(alpha_@{i}*Y_@{i}/M_@{i})^(1/epsY_@{i}) = PM_@{i};
    //% Demand for imports
        MC_@{i}*(alphaV_@{i}*Y_@{i}/V_@{i})^(1/epsY_@{i}) = PV;    
    //% Labor Demand
        MC_@{i}*((1-alpha_@{i}-alphaV_@{i})*Y_@{i}/L_@{i})^(1/epsY_@{i}) = PL_@{i};

    //% Market Clearing in Each Sector
        Y_@{i} = CHs_@{i} + CHg_@{i} + chiX_@{i}*X*PX/PH_@{i}
        @#for j in 1:nsec
            + beta_@{j}_@{i}*(PM_@{j}/PH_@{i})^epsM_@{j}*M_@{j}
        @#endfor
        ;

    //% Labor costs FOC  (allows for firing costs)

        PL_@{i} = w 
        + (cl_@{i}*w*L_@{i}*(L_@{i}/L_@{i}(-1)-1)/L_@{i}(-1) 
        + cl_@{i}/2*w*(L_@{i}/L_@{i}(-1)-1)^2) 
        - beta*((C(+1)/C)^-gamma)*cl_@{i}*w(+1)*(L_@{i}(+1)/L_@{i}-1)*L_@{i}(+1)^2/(L_@{i}^2);

    //% Rotemberg Pricing FOC
    1 - epsilon + epsilon*MC_@{i}/PH_@{i}
    - kappa_@{i}*(pi*PH_@{i}/PH_@{i}(-1)-1)*pi*PH_@{i}/PH_@{i}(-1) 
    + beta*(C(+1)/C)^(-gamma)*kappa_@{i}*(pi(+1)*PH_@{i}(+1)/PH_@{i}-1)
    *(pi(+1)*PH_@{i}(+1)/PH_@{i})^2/pi(+1)*Y_@{i}(+1)/Y_@{i} =0; 
    
    @#endfor


//% Taylor Rule
r = (1-rhoirule)*(1/beta) + rhoirule*r(-1) 
+ (1-rhoirule)*phi*(pi-1) + vi;

//% Shock processes
vi = rhoi*vi(-1) + sigma_i*eps_i;

//% Goods/Services Demand Process
exp(om_g) = (1+rho_om2-rho_om1)*ombar + rho_om1*exp(om_g(-1)) - rho_om2*exp(om_g(-2)) + sigma_om*eps_om;
exp(om_s) = (1+rho_om2-rho_om1)*(1-ombar) + rho_om1*exp(om_s(-1)) - rho_om2*exp(om_s(-2))- sigma_om*eps_om;

chi = (1-rho) + rho*chi(-1) + sigma_L_agg*epschi;

@#for i in 1:nsec
    //% TFP
    exp(A_@{i}) = (1+rho_tfp2-rho_tfp1) + rho_tfp1*exp(A_@{i}(-1)) - rho_tfp2*exp(A_@{i}(-2)) + isigma_tfp_@{i}*epsA_@{i};
@#endfor

//% Small open economy equations
C^(-gamma) = beta*(C(+1)^-gamma)*r_star*pi_e(+1)/pi(+1);

r_star = Rworld*exp(-chii_b*(bbar-Q*Bstar/GDP));

Rworld = Pistar_ss/beta;  // Rworld must satisfy Euler equation in steady state

pi_e                 =  Q/Q(-1)*pi/Pipstar; 


Q*Bstar =-TB+r_star(-1)*Q*Bstar(-1)/Pipstar;

PX=(1
    @#for i in 1:nsec
       *PH_@{i}^chiX_@{i}
    @#endfor
    );

TB = PX*X-PV*IMP;

X = omegaX*(PX/Q)^(-etastar)*Ystar;

Ystar = ystar_ss*exp(xi_psi);
xi_psi = rho_psi*xi_psi(-1) + sigma_psi*eps_psi;

%1 - epsilonV + epsilonV*Q*PVstar/PV 
%- kappaV*(pi*PV/PV(-1)-1)*pi*PV/PV(-1)
%    + beta*(C(+1)/C)^(-gamma)*kappaV*(pi(+1)*PV(+1)/PV-1)
%    *(pi(+1)*PV(+1)/PV)^2/pi(+1)*V(+1)/V =0; 
    

PV=Q*PVstar;


mkupV=1; //PV/(Q*PVstar);

IMP=mkupV*(V+CFs+CFg);

V=(0
    @#for i in 1:nsec
        +V_@{i}
    @#endfor
    );

CFs=(0
    @#for i in 1:nsec
        + CFs_@{i}
    @#endfor
    );

CFg=(0
    @#for i in 1:nsec
        + CFg_@{i}
    @#endfor
    );

Pipstar= Pistar_ss;
PVstar=PVstar_ss;

end;

initval;
eps_om  = 0;
xi_psi  = 0;
pi_g    = 1;
pi_s = 1;
pi = pi_ss;
chi = 1;
r = r_ss;
om_g=ombar;
om_s=1-ombar;

r_star=Rworld_ss; 
Rworld=Rworld_ss;  
pi_e =1;
Ystar=ystar_ss;
Bstar = Bstar_ss;  // Use computed steady state debt, not bbar!
Pipstar=Pistar_ss;
PVstar=PVstar_ss;
C = C_ss;
C_g =C_g_ss;
C_s =C_s_ss;

@#for i in 1:nsec
    Cg_@{i}  = Cgi_ss@{i};
    Cs_@{i}  = Csi_ss@{i};
    CFg_@{i} = CFg_ss@{i};
    CFs_@{i} = CFs_ss@{i};
    CHg_@{i} = CHg_ss@{i};
    CHs_@{i} = CHs_ss@{i};
    V_@{i}   = Vi_ss@{i};
    PH_@{i}  = PH_ss@{i};
    MC_@{i}  = MC_ss@{i};
    Y_@{i}   = Y_ss@{i};
    L_@{i}   = L_ss@{i};
    P_@{i}   = P_ss@{i};
    PM_@{i}  = PMi_ss@{i};
    M_@{i}   = Mi_ss@{i};
    A_@{i}   = 0;

@#endfor
N = N_ss;
p_s = p_s_ss;
p_g = p_g_ss;
mkupV = epsilonV/(epsilonV-1);
Q    = Q_ss;
PX   = PX_ss;
X = omegaX*(PX/Q)^(-etastar)*Ystar;
V = V_ss;
PV = Q*PVstar*mkupV;
TB = PX*X-PV*(CF_ss+V_ss);
w  = w_ss;
GDP = C+TB;
@#for i in 1:nsec
    PL_@{i} = PL_ss@{i};
@#endfor
Lab_costs = 0;
IMP = IMP_ss;
end;

steady;
ss1 = oo_.steady_state;


% Compute and display residuals
fprintf('\n--- Computing Residuals of Static Equations ---\n');
resid;


% STOP HERE TO REVIEW RESIDUALS
%error('STOPPING: Review residuals above. Comment this line to continue.');



// Shocks
shocks;
var eps_om=1.0;
var eps_i=1.0;
var epschi=0.0;
var eps_psi=1.0;
@#for z in 1:nsec
   var epsA_@{z}=1.0;
@#endfor
end;

check;

// Run stochastic simulation (add periods to produce oo_.endo_simul)
stoch_simul(order=1, irf=150, periods=200, nograph);

//@#include "solution_block.mod" 