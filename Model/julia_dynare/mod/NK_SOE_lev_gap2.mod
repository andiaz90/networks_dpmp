//% N sector NK model with Input-Output Linkages

// Include blocks defined elsewhere
@#include "definition_block_nsec.mod" 
@#include "definition_block_io.mod" 
@#include "definition_block_lab.mod" 

// Define variables
var chi xi Lab_costs Price_costs w VA C Ctot Ctotg Ctots M C_g C_s p_g p_s N pi pi_g pi_s r om_g om_s vi Y
PV r_star Rworld pi_e Bstar Pipstar Q TB PX X Ystar V PVstar GDP IMP mkupV CFs CFg
POstar PO VOil VNon

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
    @#for i in 1:nsec
    PIV_@{i}
    VOil_@{i}
    VNon_@{i}
    @#endfor
    ;

var C_f C_g_f C_s_f CFs_f CFg_f PV_f p_g_f p_s_f r_f N_f Lab_costs_f
w_f Y_f  Ctot_f Ctots_f Ctotg_f VA_f GDP_f TB_f M_f r_star_f pi_e_f Q_f Bstar_f
PX_f X_f IMP_f V_f Y_g_f Y_s_f PO_f VOil_f VNon_f

@#for j in 1:nsec
    P_f_@{j}
    PM_f_@{j}
    PH_f_@{j}
    Cg_f_@{j}
    Cs_f_@{j}
    CHg_f_@{j}
    CFg_f_@{j}
    CHs_f_@{j}
    CFs_f_@{j}
    L_f_@{j}
    Y_f_@{j}
    M_f_@{j}
    V_f_@{j}
    MC_f_@{j}
    PL_f_@{j}
    PIV_f_@{j}
    VOil_f_@{j}
    VNon_f_@{j}
  @#endfor
  ;

// Gap variables (log-deviations of NK from flex-price equilibrium)
var Ygap Ngap GDPgap
    Cgap_g Cgap_s
    Ygap_g Ygap_s
    Y_g Y_s
    @#for i in 1:nsec
        Ygap_@{i}
        Lgap_@{i}
    @#endfor
    ;


varexo eps_om eps_i epschi eps_pvstar eps_postar
    @#for i in 1:nsec
        epsA_@{i}
    @#endfor
    eps_xi
;


parameters gamma psi beta phi epsilon rho rho_om1 rho_tfp1 rho_om2 rho_tfp2 ombar
rhoi rhoirule ilabcosts sigma_i sigma_om sigma_L_agg Rworld_ss
bbar chii_b omegaX epsilonX ystar_ss etastar kappaV epsilonV sigmaH Pistar_ss PVstar_ss
rho_pvstar sigma_pvstar rho_xi sigma_xi
// Oil sector parameters
epsilonV_oil rho_postar sigma_postar POstar_ss shock_eps_postar
Ctot_ss Ctotg_ss Ctots_ss VA_ss M_tot_ss Y_ss IMP_ss
// Shock activation parameters (set by params_jl.mod; 0=off, 1=on)
shock_eps_om shock_eps_i shock_eps_pvstar shock_eps_xi
// Scalar SS values used in initval block (set by params_jl.mod)
pi_ss r_ss w_ss N_ss GDP_ss C_ss C_g_ss C_s_ss p_g_ss p_s_ss
Bstar_ss Q_ss TB_ss PX_ss V_ss CF_ss CFg_total_ss CFs_total_ss

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
        shock_epsA_@{i}
        alphaOilShare_@{i}
        PIV_ss@{i}
        // Sectoral SS values used in initval block
        PH_ss@{i}
        MC_ss@{i}
        Y_ss@{i}
        L_ss@{i}
        Cgi_ss@{i}
        Csi_ss@{i}
        P_ss@{i}
        PMi_ss@{i}
        Mi_ss@{i}
        CFg_ss@{i}
        CFs_ss@{i}
        CHg_ss@{i}
        CHs_ss@{i}
        Vi_ss@{i}
    @#endfor
    ;

      
// -----------------------------------------------------------------------
// All parameter values are loaded from a Julia-generated include file.
// main_SOE_gap.jl writes mod/params_jl.mod before calling @dynare.
// -----------------------------------------------------------------------
@#include "params_jl.mod"



model(bytecode);

//% Overall Consumption Basket

//% Goods/services aggregator
C = (C_g/exp(om_g))^exp(om_g)*(C_s/exp(om_s))^exp(om_s) ;

//% Goods/services demands
p_g = exp(om_g)*C/C_g;
p_s = exp(om_s)*C/C_s;



//% Goods/services aggregator
C_f = (C_g_f/exp(om_g))^exp(om_g)*(C_s_f/exp(om_s))^exp(om_s) ;

//% Goods/services demands
p_g_f = exp(om_g)*C_f/C_g_f;
p_s_f = exp(om_s)*C_f/C_s_f;



% Good j price index.

@#for j in 1:nsec
    P_@{j}=(varrho_@{j}^(sigmaH)*PH_@{j}^(1-sigmaH)+(1-varrho_@{j})^(sigmaH)*PV^(1-sigmaH))^(1/(1-sigmaH));
@#endfor

@#for j in 1:nsec
    P_f_@{j}=(varrho_@{j}^(sigmaH)*PH_f_@{j}^(1-sigmaH)+(1-varrho_@{j})^(sigmaH)*PV_f^(1-sigmaH))^(1/(1-sigmaH));
@#endfor


//% Relative Price of Goods
p_g = (1
@#for j in 1:nsec
    *(P_@{j})^gammag_@{j}
@#endfor
);

p_g_f = (1
@#for j in 1:nsec
    *(P_f_@{j})^gammag_@{j}
@#endfor
);


//% Relative Price of Services
p_s = (1
@#for j in 1:nsec
    *(P_@{j})^gammas_@{j}
@#endfor
);

//% Relative Price of Services
p_s_f = (1
@#for j in 1:nsec
    *(P_f_@{j})^gammas_@{j}
@#endfor
);

@#for j in 1:nsec
    Cg_@{j}=gammag_@{j}*(p_g/P_@{j})*C_g;
    Cs_@{j}=gammas_@{j}*(p_s/P_@{j})*C_s;
@#endfor

@#for j in 1:nsec
    Cg_f_@{j}=gammag_@{j}*(p_g_f/P_f_@{j})*C_g_f;
    Cs_f_@{j}=gammas_@{j}*(p_s_f/P_f_@{j})*C_s_f;
@#endfor



% Demands for goods
@#for j in 1:nsec   
    CHg_@{j} =    varrho_@{j}^(sigmaH)*(PH_@{j}/P_@{j})^(-sigmaH)*Cg_@{j};
    CFg_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV/P_@{j})^(-sigmaH)*Cg_@{j};    
    CHs_@{j} =    varrho_@{j}^(sigmaH)*(PH_@{j}/P_@{j})^(-sigmaH)*Cs_@{j};
    CFs_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV/P_@{j})^(-sigmaH)*Cs_@{j};
@#endfor

% Demands for goods
@#for j in 1:nsec   
    CHg_f_@{j} =    varrho_@{j}^(sigmaH)*(PH_f_@{j}/P_f_@{j})^(-sigmaH)*Cg_f_@{j};
    CFg_f_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV_f/P_f_@{j})^(-sigmaH)*Cg_f_@{j};    
    CHs_f_@{j} =    varrho_@{j}^(sigmaH)*(PH_f_@{j}/P_f_@{j})^(-sigmaH)*Cs_f_@{j};
    CFs_f_@{j}      = (1-varrho_@{j})^(sigmaH)*(PV_f/P_f_@{j})^(-sigmaH)*Cs_f_@{j};
@#endfor



//% Household Euler Equation (with preference/demand shock xi)
//% xi > 1: more impatient, want to consume now -> demand shock -> positive GDP-pi comovement
xi*C^(-gamma) = beta*xi(+1)*(C(+1)^-gamma)*r/pi(+1);

//% Household Euler Equation flex-price (same preference shock applies)
xi*C_f^(-gamma) = beta*xi(+1)*(C_f(+1)^-gamma)*r_f;

//% Relative Inflation: Goods vs Services
pi_g*C_g/C_g(-1)*(exp(om_s)/exp(om_s(-1))) = pi_s*C_s/C_s(-1)*(exp(om_g)/exp(om_g(-1)));

//% Labor Market Clearing
N = (
    @#for i in 1:nsec
        +L_@{i} 
        +cl_@{i}/2*L_@{i}*(L_@{i}/L_@{i}(-1)-1)^2
    @#endfor
    );
//% Labor Market Clearing
N_f = (
    @#for i in 1:nsec
        +L_f_@{i} 
        +cl_@{i}/2*L_f_@{i}*(L_f_@{i}/L_f_@{i}(-1)-1)^2
    @#endfor
    );


Lab_costs = (
    @#for i in 1:nsec
        + ilabcosts*cl_@{i}/2*L_@{i}*(L_@{i}/L_@{i}(-1)-1)^2
    @#endfor
    );

Lab_costs_f = (
    @#for i in 1:nsec
        + ilabcosts*cl_@{i}/2*L_f_@{i}*(L_f_@{i}/L_f_@{i}(-1)-1)^2
    @#endfor
    );
    
//% Labor leisure condition
C^(-gamma)*w = chi*N^psi;

//% Labor leisure condition
C_f^(-gamma)*w_f = chi*N_f^psi;


//% Relative Inflation: Goods vs Overall
pi_g*C_g/C_g(-1) = exp(om_g)/exp(om_g(-1))*pi*C/C(-1);


//% Gross Output
Y = (
    @#for i in 1:nsec
        +Y_@{i}
    @#endfor
    );

//% Gross Output
Y_f = (
    @#for i in 1:nsec
        +Y_f_@{i}
    @#endfor
    );



//% Total Consumption (Note: = Value Added)
Ctot = (
@#for i in 1:nsec
    +gammag_@{i}*C_g*p_g/P_@{i} + gammas_@{i}*C_s*p_s/P_@{i}
@#endfor
);

//% Total Consumption (Note: = Value Added)
Ctot_f = (
@#for i in 1:nsec
    +gammag_@{i}*C_g_f*p_g_f/P_f_@{i} + gammas_@{i}*C_s_f*p_s_f/P_f_@{i}
@#endfor
);


//% Total Consumption: Services
Ctots = (
@#for i in 1:nsec
    +gammas_@{i}*C_s*p_s/P_@{i}
@#endfor
);

//% Total Consumption: Services
Ctots_f = (
@#for i in 1:nsec
    +gammas_@{i}*C_s_f*p_s_f/P_f_@{i}
@#endfor
);

//% Total Consumption: Goods
Ctotg = (
@#for i in 1:nsec
    +gammag_@{i}*C_g*p_g/P_@{i}
@#endfor
);


//% Total Consumption: Goods
Ctotg_f = (
@#for i in 1:nsec
    +gammag_@{i}*C_g_f*p_g_f/P_f_@{i}
@#endfor
);

//% Value Added
VA = (
    @#for i in 1:nsec
        +Y_@{i} - M_@{i}
    @#endfor
);

//% Value Added
VA_f = (
    @#for i in 1:nsec
        +Y_f_@{i} - M_f_@{i}
    @#endfor
);


GDP = C + TB;
GDP_f = C_f + TB_f;


//% Total Intermediate Use
M = (
@#for i in 1:nsec
    @#for j in 1:nsec
        + beta_@{j}_@{i}*(PM_@{j}/PH_@{i})^epsM_@{j}*M_@{j}
    @#endfor
@#endfor
);

//% Total Intermediate Use
M_f = (
@#for i in 1:nsec
    @#for j in 1:nsec
        + beta_@{j}_@{i}*(PM_f_@{j}/PH_f_@{i})^epsM_@{j}*M_f_@{j}
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

    //% Intermediates Price Index
    PM_f_@{i} = (
    @#for j in 1:nsec
        + beta_@{i}_@{j}*(PH_f_@{j})^(1-epsM_@{i})
    @#endfor
    )^(1/(1-epsM_@{i})); 


    //% Production Function
        Y_@{i} = exp(A_@{i})*((alpha_@{i})^(1/epsY_@{i})*(M_@{i})^((epsY_@{i}-1)/epsY_@{i}) 
                               + (alphaV_@{i})^(1/epsY_@{i})*(V_@{i})^((epsY_@{i}-1)/epsY_@{i})
                               + (1-alpha_@{i}-alphaV_@{i})^(1/epsY_@{i})*(L_@{i})^((epsY_@{i}-1)/epsY_@{i}))^(epsY_@{i}/(epsY_@{i}-1));

    //% Production Function
        Y_f_@{i} = exp(A_@{i})*((alpha_@{i})^(1/epsY_@{i})*(M_f_@{i})^((epsY_@{i}-1)/epsY_@{i}) 
                               + (alphaV_@{i})^(1/epsY_@{i})*(V_f_@{i})^((epsY_@{i}-1)/epsY_@{i})
                               + (1-alpha_@{i}-alphaV_@{i})^(1/epsY_@{i})*(L_f_@{i})^((epsY_@{i}-1)/epsY_@{i}))^(epsY_@{i}/(epsY_@{i}-1));

    //% Intermediates Demand
        MC_@{i}*(alpha_@{i}*Y_@{i}/M_@{i})^(1/epsY_@{i}) = PM_@{i};
    //% Intermediates Demand
        MC_f_@{i}*(alpha_@{i}*Y_f_@{i}/M_f_@{i})^(1/epsY_@{i}) = PM_f_@{i};

    //% Demand for imports (composite intermediate imports V_@{i} at sector-specific price PIV_@{i})
        MC_@{i}*(alphaV_@{i}*Y_@{i}/V_@{i})^(1/epsY_@{i}) = PIV_@{i};
    //% Demand for imports (flex-price)
        MC_f_@{i}*(alphaV_@{i}*Y_f_@{i}/V_f_@{i})^(1/epsY_@{i}) = PIV_f_@{i};

    //% Sector-specific composite import price index (CES between oil and non-oil imports)
    //% PIV^(1-eps) = alphaOil * PO^(1-eps) + (1-alphaOil) * PV^(1-eps)
    PIV_@{i}^(1-epsilonV_oil) = alphaOilShare_@{i}*PO^(1-epsilonV_oil) + (1-alphaOilShare_@{i})*PV^(1-epsilonV_oil);
    PIV_f_@{i}^(1-epsilonV_oil) = alphaOilShare_@{i}*PO_f^(1-epsilonV_oil) + (1-alphaOilShare_@{i})*PV_f^(1-epsilonV_oil);

    //% Oil and non-oil import demands (CES factor demands within composite V_@{i})
    VOil_@{i} = alphaOilShare_@{i}*(PIV_@{i}/PO)^epsilonV_oil*V_@{i};
    VNon_@{i} = (1-alphaOilShare_@{i})*(PIV_@{i}/PV)^epsilonV_oil*V_@{i};
    VOil_f_@{i} = alphaOilShare_@{i}*(PIV_f_@{i}/PO_f)^epsilonV_oil*V_f_@{i};
    VNon_f_@{i} = (1-alphaOilShare_@{i})*(PIV_f_@{i}/PV_f)^epsilonV_oil*V_f_@{i};

    //% Labor Demand
        MC_@{i}*((1-alpha_@{i}-alphaV_@{i})*Y_@{i}/L_@{i})^(1/epsY_@{i}) = PL_@{i};
    //% Labor Demand
        MC_f_@{i}*((1-alpha_@{i}-alphaV_@{i})*Y_f_@{i}/L_f_@{i})^(1/epsY_@{i}) = PL_f_@{i};


    //% Market Clearing in Each Sector
        Y_@{i} = CHs_@{i} + CHg_@{i} + chiX_@{i}*X*PX/PH_@{i}
        @#for j in 1:nsec
            + beta_@{j}_@{i}*(PM_@{j}/PH_@{i})^epsM_@{j}*M_@{j}
        @#endfor
        ;
    //% Market Clearing in Each Sector
        Y_f_@{i} = CHs_f_@{i} + CHg_f_@{i} + chiX_@{i}*X_f*PX_f/PH_f_@{i}
        @#for j in 1:nsec
            + beta_@{j}_@{i}*(PM_f_@{j}/PH_f_@{i})^epsM_@{j}*M_f_@{j}
        @#endfor
        ;


    //% Labor costs FOC  (allows for firing costs)
        PL_@{i} = w 
        + (cl_@{i}*w*L_@{i}*(L_@{i}/L_@{i}(-1)-1)/L_@{i}(-1) 
        + cl_@{i}/2*w*(L_@{i}/L_@{i}(-1)-1)^2) 
        - beta*((C(+1)/C)^-gamma)*cl_@{i}*w(+1)*(L_@{i}(+1)/L_@{i}-1)*L_@{i}(+1)^2/(L_@{i}^2);

    //% Labor costs FOC  (allows for firing costs)

        PL_f_@{i} = w_f 
        + (cl_@{i}*w_f*L_f_@{i}*(L_f_@{i}/L_f_@{i}(-1)-1)/L_f_@{i}(-1) 
        + cl_@{i}/2*w_f*(L_f_@{i}/L_f_@{i}(-1)-1)^2) 
        - beta*((C_f(+1)/C_f)^-gamma)*cl_@{i}*w_f(+1)*(L_f_@{i}(+1)/L_f_@{i}-1)*L_f_@{i}(+1)^2/(L_f_@{i}^2);



    //% Rotemberg Pricing FOC
    //% Derivation: adjustment cost (kappa/2)*(Pi_H - 1)^2 * Y, where Pi_H = pi*PH/PH(-1)
    //% FOC with respect to PH_@{i} yields the standard Rotemberg NKPC:
    //% 1 - eps + eps*MC/PH - kappa*(Pi_H-1)*Pi_H + beta*M_{t+1}*kappa*(Pi_H(+1)-1)*Pi_H(+1)*Y(+1)/Y = 0
    1 - epsilon + epsilon*MC_@{i}/PH_@{i}
    - kappa_@{i}*(pi*PH_@{i}/PH_@{i}(-1)-1)*pi*PH_@{i}/PH_@{i}(-1)
    + beta*(C(+1)/C)^(-gamma)*kappa_@{i}*(pi(+1)*PH_@{i}(+1)/PH_@{i}-1)
    *(pi(+1)*PH_@{i}(+1)/PH_@{i})*Y_@{i}(+1)/Y_@{i} =0;

    //% Rotemberg Pricing FOC
    1 - epsilon + epsilon*MC_f_@{i}/PH_f_@{i}=0; 



    @#endfor


//% Aggregate oil and non-oil intermediate imports
VOil = (0
@#for i in 1:nsec
    + VOil_@{i}
@#endfor
);
VNon = (0
@#for i in 1:nsec
    + VNon_@{i}
@#endfor
);
VOil_f = (0
@#for i in 1:nsec
    + VOil_f_@{i}
@#endfor
);
VNon_f = (0
@#for i in 1:nsec
    + VNon_f_@{i}
@#endfor
);

//% Taylor Rule
r = (1-rhoirule)*(1/beta) + rhoirule*r(-1)
+ (1-rhoirule)*phi*(pi-1) + vi;

//% Shock processes
vi = rhoi*vi(-1) + sigma_i*eps_i;

//% Goods/Services Demand Process
exp(om_g) = (1+rho_om2-rho_om1)*ombar + rho_om1*exp(om_g(-1)) - rho_om2*exp(om_g(-2)) + sigma_om*eps_om;
exp(om_s) = (1+rho_om2-rho_om1)*(1-ombar) + rho_om1*exp(om_s(-1)) - rho_om2*exp(om_s(-2))- sigma_om*eps_om;

chi = (1-rho) + rho*chi(-1) + sigma_L_agg*epschi;

//% Preference/demand shock (AR(1), SS = 1)
xi = (1-rho_xi) + rho_xi*xi(-1) + sigma_xi*eps_xi;

@#for i in 1:nsec
    //% TFP
    exp(A_@{i}) = (1+rho_tfp2-rho_tfp1) + rho_tfp1*exp(A_@{i}(-1)) - rho_tfp2*exp(A_@{i}(-2)) + isigma_tfp_@{i}*epsA_@{i};
@#endfor

//% Small open economy equations
C^(-gamma) = beta*(C(+1)^-gamma)*r_star*pi_e(+1)/pi(+1);

//% Small open economy equations
C_f^(-gamma) = beta*(C_f(+1)^-gamma)*r_star_f*pi_e_f(+1);


r_star = Rworld*exp(-chii_b*(bbar-Q*Bstar/GDP));

r_star_f = Rworld*exp(-chii_b*(bbar-Q_f*Bstar_f/GDP_f));

Rworld = Pistar_ss/beta;  // Rworld must satisfy Euler equation in steady state

pi_e                 =  Q/Q(-1)*pi/Pipstar; 

pi_e_f                 =  Q_f/Q_f(-1); 


Q*Bstar =-TB+r_star(-1)*Q*Bstar(-1)/Pipstar;

Q_f*Bstar_f =-TB_f+r_star_f(-1)*Q_f*Bstar_f(-1);


PX=(1
    @#for i in 1:nsec
       *PH_@{i}^chiX_@{i}
    @#endfor
    );

PX_f=(1
    @#for i in 1:nsec
       *PH_f_@{i}^chiX_@{i}
    @#endfor
    );

//% Trade balance: exports minus oil imports (at PO) and non-oil imports (at PV)
//% Consumer imports (CFs, CFg) priced at PV (non-oil)
TB = PX*X - PO*VOil - PV*(VNon+CFs+CFg);
TB_f = PX_f*X_f - PO_f*VOil_f - PV_f*(VNon_f+CFs_f+CFg_f);

X = omegaX*(PX/Q)^(-etastar)*Ystar;
X_f = omegaX*(PX_f/Q_f)^(-etastar)*Ystar;

Ystar=ystar_ss;

%1 - epsilonV + epsilonV*Q*PVstar/PV 
%- kappaV*(pi*PV/PV(-1)-1)*pi*PV/PV(-1)
%    + beta*(C(+1)/C)^(-gamma)*kappaV*(pi(+1)*PV(+1)/PV-1)
%    *(pi(+1)*PV(+1)/PV)^2/pi(+1)*V(+1)/V =0; 
    

PV=Q*PVstar;
PV_f=Q_f*PVstar;


mkupV=1; //PV/(Q*PVstar);

IMP=mkupV*(V+CFs+CFg);
IMP_f=mkupV*(V_f+CFs_f+CFg_f);

V=(0
    @#for i in 1:nsec
        +V_@{i}
    @#endfor
    );

V_f=(0
    @#for i in 1:nsec
        +V_f_@{i}
    @#endfor
    );    

CFs=(0
    @#for i in 1:nsec
        + CFs_@{i}
    @#endfor
    );

CFs_f=(0
    @#for i in 1:nsec
        + CFs_f_@{i}
    @#endfor
    );

CFg=(0
    @#for i in 1:nsec
        + CFg_@{i}
    @#endfor
    );

CFg_f=(0
    @#for i in 1:nsec
        + CFg_f_@{i}
    @#endfor
    );

Pipstar= Pistar_ss;
log(PVstar/PVstar_ss) = rho_pvstar*log(PVstar(-1)/PVstar_ss) + sigma_pvstar*eps_pvstar;

//% Oil price process (world oil price in foreign currency, pass-through via exchange rate)
PO = Q*POstar;
PO_f = Q_f*POstar;
log(POstar/POstar_ss) = rho_postar*log(POstar(-1)/POstar_ss) + sigma_postar*eps_postar;

//% Output and Employment Gaps (log-deviations from flex-price equilibrium)
Ygap    = log(Y)   - log(Y_f);
Ngap    = log(N)   - log(N_f);
GDPgap  = log(GDP) - log(GDP_f);

//% Sectoral output aggregates: goods and services (NK)
Y_g = (0
@#for i in 1:nsec
    + dummyg_@{i}*Y_@{i}
@#endfor
);

Y_s = (0
@#for i in 1:nsec
    + dummys_@{i}*Y_@{i}
@#endfor
);

//% Sectoral output aggregates: goods and services (flex-price)
Y_g_f = (0
@#for i in 1:nsec
    + dummyg_@{i}*Y_f_@{i}
@#endfor
);

Y_s_f = (0
@#for i in 1:nsec
    + dummys_@{i}*Y_f_@{i}
@#endfor
);

//% Proper Output Gaps: goods and services sectors
Ygap_g  = log(Y_g)   - log(Y_g_f);
Ygap_s  = log(Y_s)   - log(Y_s_f);

//% Consumption Gaps: goods and services (log-deviation of NK from flex-price)
Cgap_g  = log(Ctotg) - log(Ctotg_f);
Cgap_s  = log(Ctots) - log(Ctots_f);

@#for i in 1:nsec
    Ygap_@{i} = log(Y_@{i})   - log(Y_f_@{i});
    Lgap_@{i} = log(L_@{i})   - log(L_f_@{i});
@#endfor

end;

initval;
eps_om = 0;
eps_pvstar = 0;
eps_xi = 0;
pi_g = 1;
pi_s = 1;
pi = pi_ss;
chi = 1;
xi = 1;
r = r_ss;
r_f = r_ss;
vi = 0;

om_g=log(ombar);     // SS: exp(om_g) = ombar, so om_g = log(ombar)
om_s=log(1-ombar);  // SS: exp(om_s) = 1-ombar

r_star=Rworld_ss; 
r_star_f=Rworld_ss; 

Rworld=Rworld_ss;  
pi_e =1;
pi_e_f =1;

Ystar=ystar_ss;
Bstar = Bstar_ss;
Bstar_f = Bstar_ss;

Pipstar=Pistar_ss;
PVstar=PVstar_ss;
C = C_ss;
C_g =C_g_ss;
C_s =C_s_ss;
C_f = C_ss;
C_g_f = C_g_ss;
C_s_f = C_s_ss;
p_g_f = p_g_ss;
p_s_f = p_s_ss;
Price_costs = 0;
Lab_costs = 0;

@#for i in 1:nsec
    Cg_@{i}  = Cgi_ss@{i};
    Cs_@{i}  = Csi_ss@{i};    
    CFg_@{i} = CFg_ss@{i};
    CFs_@{i} = CFs_ss@{i};

    Cg_f_@{i}  = Cgi_ss@{i};
    Cs_f_@{i}  = Csi_ss@{i};    
    CFg_f_@{i} = CFg_ss@{i};
    CFs_f_@{i} = CFs_ss@{i};

    CHg_@{i} = CHg_ss@{i};
    CHs_@{i} = CHs_ss@{i};
    V_@{i}   = Vi_ss@{i};
    PH_@{i}  = PH_ss@{i};

    CHg_f_@{i} = CHg_ss@{i};
    CHs_f_@{i} = CHs_ss@{i};
    V_f_@{i}   = Vi_ss@{i};
    PH_f_@{i}  = PH_ss@{i};

    MC_@{i}  = MC_ss@{i};
    Y_@{i}   = Y_ss@{i};
    L_@{i}   = L_ss@{i};
    P_@{i}   = P_ss@{i};

    MC_f_@{i}  = MC_ss@{i};
    Y_f_@{i}   = Y_ss@{i};
    L_f_@{i}   = L_ss@{i};
    P_f_@{i}   = P_ss@{i};

    PM_@{i}  = PMi_ss@{i};
    M_@{i}   = Mi_ss@{i};
    PL_@{i} = PL_ss@{i};    

    PM_f_@{i}  = PMi_ss@{i};
    M_f_@{i}   = Mi_ss@{i};
    PL_f_@{i} = PL_ss@{i};    

    A_@{i}   = 0;
@#endfor



C_f = C_ss;
C_g_f = C_g_ss;
C_s_f = C_s_ss;
N = N_ss;
N_f = N_ss;

p_s = p_s_ss;
p_g = p_g_ss;

p_s_f = p_s_ss;
p_g_f = p_g_ss;

mkupV = 1;
Q    = Q_ss;
Q_f    = Q_ss;

PX   = PX_ss;
PX_f   = PX_ss;

X = omegaX*(PX/Q)^(-etastar)*Ystar;
X_f = omegaX*(PX_f/Q_f)^(-etastar)*Ystar;

V = V_ss;
V_f = V_ss;

PV = Q*PVstar*mkupV;
PV_f = Q_f*PVstar*mkupV;

POstar = POstar_ss;
PO = Q_ss*POstar_ss;
PO_f = Q_ss*POstar_ss;

//% Oil/non-oil split at SS: since POstar_ss = PVstar_ss = 1, PO_ss = PV_ss,
//% so PIV_i_ss = PV_ss for all i (independent of alphaOilShare_i).
//% Initial values for VOil/VNon are starting guesses; Dynare's steady; will compute exact SS.
//% TB_ss = PX*X - PO_ss*VOil_ss - PV_ss*VNon_ss - PV_ss*CF_ss
//%       = PX*X - PV_ss*(VOil_ss+VNon_ss+CF_ss) = PX*X - PV*(V_ss+CF_ss)  [since PO_ss=PV_ss]
VOil = 0;
VNon = V_ss;
VOil_f = 0;
VNon_f = V_ss;

@#for i in 1:nsec
    PIV_@{i}   = PIV_ss@{i};
    PIV_f_@{i} = PIV_ss@{i};
    VOil_@{i}  = 0;
    VNon_@{i}  = Vi_ss@{i};
    VOil_f_@{i} = 0;
    VNon_f_@{i} = Vi_ss@{i};
@#endfor

//% TB_ss value is PX*X - PV*(V_ss+CF_ss) because PO_ss = PV_ss (see above)
TB = PX*X - PV*(V_ss+CF_ss);
TB_f = PX_f*X_f - PV_f*(V_ss+CF_ss);

w  = w_ss;
w_f  = w_ss;

GDP = C+TB;
GDP_f = C_f+TB_f;

Lab_costs_f = 0;

IMP = IMP_ss;
IMP_f = IMP_ss;

Y_f = Y_ss;
Y   = Y_ss;

Ctot  = Ctot_ss;
Ctotg = Ctotg_ss;
Ctots = Ctots_ss;
VA    = VA_ss;
M     = M_tot_ss;

Ctot_f  = Ctot_ss;
Ctotg_f = Ctotg_ss;
Ctots_f = Ctots_ss;
VA_f    = VA_ss;
M_f     = M_tot_ss;
CFs = CFs_total_ss;
CFg = CFg_total_ss;
CFs_f = CFs_total_ss;
CFg_f = CFg_total_ss;

Ygap   = 0;
Ngap   = 0;
GDPgap = 0;
Ygap_g = 0;
Ygap_s = 0;
Cgap_g = 0;
Cgap_s = 0;
Y_g = (0
@#for i in 1:nsec
    + dummyg_@{i}*Y_ss@{i}
@#endfor
);
Y_s = (0
@#for i in 1:nsec
    + dummys_@{i}*Y_ss@{i}
@#endfor
);
Y_g_f = (0
@#for i in 1:nsec
    + dummyg_@{i}*Y_ss@{i}
@#endfor
);
Y_s_f = (0
@#for i in 1:nsec
    + dummys_@{i}*Y_ss@{i}
@#endfor
);
@#for i in 1:nsec
    Ygap_@{i} = 0;
    Lgap_@{i} = 0;
@#endfor

end;

steady;
// ss1 = oo_.steady_state;  ← MATLAB-only: removed for Dynare.jl compatibility
// resid;                    ← MATLAB-only: removed for Dynare.jl compatibility

// Shocks
shocks;
var eps_om=shock_eps_om;
var eps_i=shock_eps_i;
var epschi=0.0;
var eps_pvstar=shock_eps_pvstar;
var eps_postar=shock_eps_postar;
var eps_xi=shock_eps_xi;
@#for z in 1:nsec
   var epsA_@{z}=shock_epsA_@{z};
@#endfor
end;

check;

// Run stochastic simulation (add periods to produce oo_.endo_simul)
stoch_simul(order=1, irf=150, periods=200, nograph);

//@#include "solution_block.mod" 