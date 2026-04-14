function loss = loss_model(x,param_name,oo_,M_,options_,nsec_val,names,p_d,y_d,l_d,goods,services,share)

load numit numit
format compact

%% ASSIGN ESTIMATION
k=1;
for i=1:numel(x)

    % For one parameter equal to one x
    if iscell(param_name{i})==0
    eval(([ param_name{i} '=x(i)' ]))
    param_name_expand{k}=param_name{i};
    k=k+1;
    end
    
    % For more than one parameter equal to the same x
    if iscell(param_name{i})==1
        for j=1:numel(param_name{i})
        eval(strjoin([ param_name{i}(j) '=x(i)' ]))
        param_name_expand{k}=strjoin(param_name{i}(j));
        k=k+1;        
        end
    end
    
end



    
%% write parameters that enter estimation into Dynare
if strmatch('sigma_L_agg_val',param_name_expand)
    M_.params(strmatch('sigma_L_agg',M_.param_names,'exact'))=sigma_L_agg_val;
end

for i=1:nsec_val
    if strmatch('modcl',param_name_expand)
    M_.params(strmatch(['cl_',num2str(i)],M_.param_names,'exact'))=modcl;
    end
    if strmatch('modclneg',param_name_expand)
    M_.params(strmatch(['clneg_',num2str(i)],M_.param_names,'exact'))=modclneg;
    end
    if strmatch('modcm',param_name_expand)
    M_.params(strmatch(['cm_',num2str(i)],M_.param_names,'exact'))=modcm;
    end
    if strmatch('modpsil',param_name_expand)
    M_.params(strmatch(['psil_',num2str(i)],M_.param_names,'exact'))=modpsil;
    end
    if strmatch('modpsim',param_name_expand)
    M_.params(strmatch(['psim_',num2str(i)],M_.param_names,'exact'))=modpsim;
    end
    if strmatch('modepsM',param_name_expand)
    M_.params(strmatch(['epsM_',num2str(i)],M_.param_names,'exact'))=modepsM;
    end
    if strmatch('modepsY',param_name_expand)
    M_.params(strmatch(['epsY_',num2str(i)],M_.param_names,'exact'))=modepsY;
    end
end




%% Call solver
[oo_.steady_state,M_.params,info] = evaluate_steady_state(oo_.steady_state,M_,options_,oo_,1);
oo_.endo_simul(:,1)=oo_.steady_state;
oo_ = perfect_foresight_solver_core(M_,options_,oo_);
    


%% Read Dynare's output

p = {};
Y = {};
A = {};
w = {};
pl = {};
L = {};
M = {};
p_m = NaN(nsec_val,1);
y_m = NaN(nsec_val,1);
l_m = NaN(nsec_val,1);
m_m = NaN(nsec_val,1);
A_m = NaN(nsec_val,1);
l_ss = NaN(nsec_val,1);

N = oo_.endo_simul(strmatch('N',M_.endo_names,'exact'),:);
r_m = oo_.endo_simul(strmatch(['r_' num2str(i)],M_.endo_names,'exact'),:)*400;

pi_g = oo_.endo_simul(strmatch('pi_g',M_.endo_names,'exact'),:);
pi_s = oo_.endo_simul(strmatch('pi_s',M_.endo_names,'exact'),:);
L_g = oo_.endo_simul(strmatch('L_g',M_.endo_names,'exact'),:);
L_s = oo_.endo_simul(strmatch('L_s',M_.endo_names,'exact'),:);
C_g = oo_.endo_simul(strmatch('C_g',M_.endo_names,'exact'),:);
C_s = oo_.endo_simul(strmatch('C_s',M_.endo_names,'exact'),:);

Ymat = [];
name_vec = strings(nsec_val,1);


for i = 1:nsec_val
    %w{i} = oo_.endo_simul(strmatch(['w_' num2str(i)],M_.endo_names,'exact'),:);
    pl{i} = oo_.endo_simul(strmatch(['PL_' num2str(i)],M_.endo_names,'exact'),:);
    p{i} = oo_.endo_simul(strmatch(['P_' num2str(i)],M_.endo_names,'exact'),:);
    L{i} = oo_.endo_simul(strmatch(['L_' num2str(i)],M_.endo_names,'exact'),:);
    Y{i} = oo_.endo_simul(strmatch(['Y_' num2str(i)],M_.endo_names,'exact'),:);
    M{i} = oo_.endo_simul(strmatch(['M_' num2str(i)],M_.endo_names,'exact'),:);
    A{i} = oo_.endo_simul(strmatch(['A_' num2str(i)],M_.endo_names,'exact'),:);
    Ymat(:,i) = Y{i};
    name_vec(i) = names{i}(1:min(50,numel(names{i})));
    A_m(i,1) = 100*(mean(A{i}(5)) - mean(A{i}(1)));
    p_m(i,1) = 100*(mean(p{i}(5)) - mean(p{i}(1)));
    y_m(i,1) = 100*(mean(Y{i}(5)) - mean(Y{i}(1)));
    l_m(i,1) = 100*(mean(L{i}(5)) - mean(L{i}(1)));
    m_m(i,1) = 100*(mean(M{i}(5)) - mean(M{i}(1)));
    l_ss(i) = exp(L{i}(1));
end




%% Calculate Moments from model and compare with data

cp=corr(p_d,p_m);
cy=corr(y_d,y_m);
cl=corr(l_d,l_m);

stdlg_d = std(l_d(goods),share(goods));
stdls_d = std(l_d(services),share(services));
stdpg_d = std(p_d(goods),share(goods));
stdps_d = std(p_d(services),share(services));
stdyg_d = std(y_d(goods),share(goods));
stdys_d = std(y_d(services),share(services));

stdlg_m = std(l_m(goods),share(goods));
stdls_m = std(l_m(services),share(services));
stdpg_m = std(p_m(goods),share(goods));
stdps_m = std(p_m(services),share(services));
stdyg_m = std(y_m(goods),share(goods));
stdys_m = std(y_m(services),share(services));

dpg = yoy(pi_g(1:20)*100);
dps = yoy(pi_s(1:20)*100);

dpg_model = mean(dpg(5)); 
dps_model = mean(dps(5)); 

lg_model = 100*(mean(L_g(5))-mean(L_g(1)));
ls_model = 100*(mean(L_s(5))-mean(L_s(1)));
l_model = 100*(mean(N(5))-mean(N(1)));

dpg_data = 8 ; % 2021Q4 relative to 2019 average
dps_data = 2 ; % 2021Q4 relative to 2019 average
l_data = -3.7; % 2021Q4 for employment relative to trend



%% Calculate loss

loss_unweighted = [ stdpg_m - stdpg_d 
    stdps_m - stdps_d
    stdyg_m - stdyg_d 
    stdys_m - stdys_d
    stdlg_m - stdlg_d 
    stdls_m - stdls_d
    cp - 1
    cy - 1
    cl - 1
    l_model - l_data
    (dpg_model-dpg_data)-(dps_model-dps_data)  ] ;

 WEIGHTS = [ 1
     1
     1
     1
     1
     1
     1
     1
     1
     1 
     1 ] ;



loss = loss_unweighted.*WEIGHTS ;
disp(loss_unweighted');

if oo_.deterministic_simulation.status~=1
    fprintf('Simulation failed')
    loss=Inf;
end



%% Save iteration
save(['lsq_loop_' num2str(numit) ])



%% Update loop
numit = numit + 1;
save numit numit

end
