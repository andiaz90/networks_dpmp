% Plot figs
load params_val modbeta modalpha

% Ensure 'other' sector classification exists
if exist('goods','var') && exist('services','var')
    try
        other = ~(goods | services);
    catch
        % Fallback if dimensions mismatch
        nsec_guess = exist('nsec_val','var') * nsec_val + (~exist('nsec_val','var')) * numel(goods);
        other = false(nsec_guess,1);
    end
else
    % If goods/services not in workspace, default to no 'other'
    if exist('nsec_val','var')
        other = false(nsec_val,1);
    else
        other = [];
    end
end

if exist('horizon_scatter_model')==0
    horizon_scatter_model=9;
    t_end=40;
    time=0:t_end-1;
end

% Extract steady states for level variables (from oo_.steady_state)
pi_idx = strmatch('pi', M_.endo_names, 'exact');
pi = oo_.endo_simul(pi_idx, :);
pi_ss = oo_.steady_state(pi_idx);

r_idx = strmatch('r', M_.endo_names, 'exact');
r = oo_.endo_simul(r_idx, :);
r_ss = oo_.steady_state(r_idx);
r_m = 400*(r - r_ss); % annual rate deviation

Ctot_idx = strmatch('Ctot', M_.endo_names, 'exact');
Ctot = oo_.endo_simul(Ctot_idx, :);
Ctot_ss = oo_.steady_state(Ctot_idx);

Ctotg_idx = strmatch('Ctotg', M_.endo_names, 'exact');
Ctotg = oo_.endo_simul(Ctotg_idx, :);
Ctotg_ss = oo_.steady_state(Ctotg_idx);

Ctots_idx = strmatch('Ctots', M_.endo_names, 'exact');
Ctots = oo_.endo_simul(Ctots_idx, :);
Ctots_ss = oo_.steady_state(Ctots_idx);

% Extract sectoral steady states from oo_.steady_state and compute scatter means
for i = 1:nsec_val
    % Get steady states from Dynare output
    p_idx = strmatch(sprintf('P_%d', i), M_.endo_names, 'exact');
    Y_idx = strmatch(sprintf('Y_%d', i), M_.endo_names, 'exact');
    L_idx = strmatch(sprintf('L_%d', i), M_.endo_names, 'exact');
    
    p_ss_i = oo_.steady_state(p_idx);
    Y_ss_i = oo_.steady_state(Y_idx);
    L_ss_i = oo_.steady_state(L_idx);
    
    % Compute percent deviations for scatter plots
    p_pct_i = 100*(p{i} - p_ss_i) / p_ss_i;
    Y_pct_i = 100*(Y{i} - Y_ss_i) / Y_ss_i;
    L_pct_i = 100*(L{i} - L_ss_i) / L_ss_i;
    
    p_m(i,1) = mean(p_pct_i(2:horizon_scatter_model));
    y_m(i,1) = mean(Y_pct_i(2:horizon_scatter_model));
    l_m(i,1) = mean(L_pct_i(2:horizon_scatter_model));
end

% Goods and services inflation series (if available)
pi_g = [];
pi_s = [];
idx_pi_g = strmatch('pi_g', M_.endo_names, 'exact');
idx_pi_s = strmatch('pi_s', M_.endo_names, 'exact');
if ~isempty(idx_pi_g)
    pi_g = oo_.endo_simul(idx_pi_g, :);
    pi_g_ss = oo_.steady_state(idx_pi_g);
end
if ~isempty(idx_pi_s)
    pi_s = oo_.endo_simul(idx_pi_s, :);
    pi_s_ss = oo_.steady_state(idx_pi_s);
end

Lab_costs = oo_.endo_simul(strmatch('Lab_costs',M_.endo_names,'exact'),:);
om_g = oo_.endo_simul(strmatch('om_g',M_.endo_names,'exact'),:);
chi_1 = oo_.endo_simul(strmatch('chi_1',M_.endo_names,'exact'),:);
w_1 = oo_.endo_simul(strmatch('w_1',M_.endo_names,'exact'),:);
vi = oo_.endo_simul(strmatch('vi',M_.endo_names,'exact'),:);
rl_m = [r_m(1) fwdmovavg(r_m(2:end),20) ];

% Calculate inflation as gross rate change -> percent (pi is gross inflation 1+rate)
pi_rate = 100*(pi./pi_ss - 1); % percent deviation from steady-state inflation
% Derive goods and services inflation if series exist; else fall back to total
if ~isempty(pi_g)
    pi_g_rate = 100*(pi_g./pi_g_ss - 1);
else
    pi_g_rate = pi_rate;
end
if ~isempty(pi_s)
    pi_s_rate = 100*(pi_s./pi_s_ss - 1);
else
    pi_s_rate = pi_rate;
end
dp_m = max(pi_rate(1:t_end));

% Create figure for 12-panel plot
figure('Name', 'Model Dynamics - 12 Panel View', 'NumberTitle', 'off', ...
       'Position', [100 100 1400 900], 'Visible', 'on');
drawnow;

subplot(3,4,1)
plot(time,100*(Ctotg(1:t_end)-Ctotg_ss)/Ctotg_ss,'b')
hold on
plot(time,100*(Ctots(1:t_end)-Ctots_ss)/Ctots_ss,'r')
hold on
plot(time,100*(Ctot(1:t_end)-Ctot_ss)/Ctot_ss,'k')
title('Consumption','fontsize',10)
legend('Goods','Serv.','Total')
ylabel('% deviation from SS')


subplot(3,4,2)
plot(time,yoy(pi_g_rate(1:t_end)),'color','b','LineWidth',1.5); hold on
plot(time,yoy(pi_s_rate(1:t_end)),'color','r','LineWidth',1.5); hold on
plot(time,yoy(pi_rate(1:t_end)),'k','LineWidth',2); hold on
title('Inflation (yoy)','fontsize',10)
fix_y_axis



subplot(3,4,3)
for i = 1:nsec_val
    % Convert p from levels to % deviation
    p_idx = strmatch(sprintf('P_%d', i), M_.endo_names, 'exact');
    p_ss_i = oo_.steady_state(p_idx);
    p_pct = 100*(p{i} - p_ss_i) / p_ss_i;
    
    if goods(i) == 1
        plot(time,p_pct(1:t_end),'color',[0, 0.4470,0.7410],'LineWidth',1.2)
    elseif services(i) == 1
        plot(time,p_pct(1:t_end),'color',[0.85,0.325,0.098],'LineWidth',1.2)
    elseif other(i) == 1
        plot(time,p_pct(1:t_end),'color',[0.3,0.3,0.3],'LineWidth',1.2)
    end
    text(2,p_pct(3),names(i),'fontsize',8)
    hold on
end
title('Sectoral Prices','fontsize',10)
ylabel('% deviation from SS')



subplot(3,4,9)
for i = 1:nsec_val
    % Convert L from levels to % deviation
    L_idx = strmatch(sprintf('L_%d', i), M_.endo_names, 'exact');
    L_ss_i = oo_.steady_state(L_idx);
    L_pct = 100*(L{i} - L_ss_i) / L_ss_i;
    
    if goods(i) == 1
        plot(time,L_pct(1:t_end),'color',[0, 0.4470, 0.7410],'LineWidth',1.2)
    elseif services(i) == 1
        plot(time,L_pct(1:t_end),'color',[0.85, 0.325, 0.098],'LineWidth',1.2)
    elseif other(i) == 1
        plot(time,L_pct(1:t_end),'color',[0.3, 0.3, 0.3],'LineWidth',1.2)
    end
    text(2,L_pct(3),names(i),'fontsize',8)
    hold on
end
title('Sectoral Labor','fontsize',10)
ylabel('% deviation from SS')
axis tight

subplot(3,4,8)
for i = 1:nsec_val
    % Convert M from levels to % deviation
    M_idx = strmatch(sprintf('M_%d', i), M_.endo_names, 'exact');
    M_ss_i = oo_.steady_state(M_idx);
    M_pct = 100*(M{i} - M_ss_i)/M_ss_i;
    
    if goods(i) == 1
        plot(time,M_pct(1:t_end),'color',[0, 0.4470,0.7410],'LineWidth',1.2)
    elseif services(i) == 1
        plot(time,M_pct(1:t_end),'color',[0.85,0.325,0.098],'LineWidth',1.2)
    elseif other(i) == 1
        plot(time,M_pct(1:t_end),'color',[0.3,0.3,0.3],'LineWidth',1.2)
    end
    text(2,M_pct(3),names(i),'fontsize',8)
    hold on
end
title('Materials','fontsize',10)
ylabel('% deviation from SS')
axis tight



subplot(3,4,5)
for i = 1:nsec_val
    % Convert Y from levels to % deviation
    Y_idx = strmatch(sprintf('Y_%d', i), M_.endo_names, 'exact');
    Y_ss_i = oo_.steady_state(Y_idx);
    Y_pct = 100*(Y{i} - Y_ss_i) / Y_ss_i;
    
    if goods(i) == 1
        plot(time,Y_pct(1:t_end),'color',[0, 0.4470,0.7410],'LineWidth',1.2)
    elseif services(i) == 1
        plot(time,Y_pct(1:t_end),'color',[0.85,0.325,0.098],'LineWidth',1.2)
    elseif other(i) == 1
        plot(time,Y_pct(1:t_end),'color',[0.3,0.3,0.3],'LineWidth',1.2)
    end
    text(2,Y_pct(3),names(i),'fontsize',8)
    hold on
end
title('Sectoral Output','fontsize',10)
ylabel('% deviation from SS')
axis tight


subplot(3,4,6)
% Plot aggregate labor (N)
N_idx = strmatch('N', M_.endo_names, 'exact');
if ~isempty(N_idx)
    N = oo_.endo_simul(N_idx, :);
    N_ss = oo_.steady_state(N_idx);
    N_pct = 100*(N - N_ss) / N_ss;
    plot(time, N_pct(1:t_end), 'k', 'LineWidth', 2)
    title('Aggregate Labor (N)', 'fontsize', 10)
    ylabel('% deviation from SS')
    grid on
    axis tight
else
    % Fallback: plot sum of sectoral labor
    L_total = zeros(size(time));
    for i = 1:nsec_val
        L_idx = strmatch(sprintf('L_%d', i), M_.endo_names, 'exact');
        L_total = L_total + oo_.endo_simul(L_idx, 1:t_end);
    end
    plot(time, L_total, 'k', 'LineWidth', 2)
    title('Total Labor Supply', 'fontsize', 10)
    ylabel('Level')
    grid on
    axis tight
end

subplot(3,4,7)
for i = 1:nsec_val
    % A_{i} is log(TFP), so convert to levels first: TFP = exp(A_{i})
    TFP_level = exp(A{i});
    
    % Get steady state TFP in levels
    A_idx = strmatch(sprintf('A_%d', i), M_.endo_names, 'exact');
    if ~isempty(A_idx)
        A_ss_log = oo_.steady_state(A_idx);
        TFP_ss = exp(A_ss_log);
    else
        % If not found, use first period as steady state
        TFP_ss = TFP_level(1);
    end
    
    % Convert to % deviation from steady state
    A_pct = 100*(TFP_level - TFP_ss)/TFP_ss;
    
    if goods(i) == 1
        plot(time,A_pct(1:t_end),'color',[0, 0.4470,0.7410],'LineWidth',1.2)
    elseif services(i) == 1
        plot(time,A_pct(1:t_end),'color',[0.85,0.325,0.098],'LineWidth',1.2)
    elseif other(i) == 1
        plot(time,A_pct(1:t_end),'color',[0.3,0.3,0.3],'LineWidth',1.2)
    end
    text(2,A_pct(2),names(i),'fontsize',8)
    hold on
end
title('Sectoral TFP','fontsize',10)
ylabel('% deviation from SS')
axis tight

[U,U_g,U_s,U_g_seb,leontief] = upstream(modbeta,modalpha,goods,M_,oo_);



subplot(3,4,4)
scatter(U(goods),p_m(goods),30,'filled','blue')
hold on
scatter(U(services),p_m(services),30,'filled','red')
hold on
scatter(U(other),p_m(other),30,'filled','black')
hold on
for i = 1:nsec_val
    text(U(i),p_m(i),names(i),'fontsize',8)
end
xlabel('Upstream tot')
ylabel('Model prices')


subplot(3,4,10)
scatter(p_d(goods),p_m(goods),30,'filled','blue')
hold on
scatter(p_d(services),p_m(services),30,'filled','red')
hold on
scatter(p_d(other),p_m(other),30,'filled','black')
hold on
for i = 1:nsec_val
    text(p_d(i),p_m(i),names(i),'fontsize',8)
end
xlabel('Prices: Data')
ylabel('Model')
grid on


subplot(3,4,11)
scatter(y_d(goods),y_m(goods),30,'filled','blue')
hold on
scatter(y_d(services),y_m(services),30,'filled','red')
hold on
scatter(y_d(other),y_m(other),30,'filled','black')
hold on
for i = 1:nsec_val
    text(y_d(i),y_m(i),names(i),'fontsize',8)
end
xlabel('Output: Data')
ylabel('Model')
grid on

subplot(3,4,12)
scatter(l_d(goods),l_m(goods),30,'filled','blue')
hold on
scatter(l_d(services),l_m(services),30,'filled','red')
hold on
scatter(l_d(other),l_m(other),30,'filled','black')
hold on
for i = 1:nsec_val
    text(l_d(i),l_m(i),names(i),'fontsize',8)
end
xlabel('Labor: Data')
ylabel('Model')
grid on

% Save the 12-panel figure
try
    print(gcf, '-dpng', '-r300', 'model_12panel_dynamics.png');
    fprintf('✓ 12-panel dynamics figure saved: model_12panel_dynamics.png\n');
catch ME
    warning('Failed to save 12-panel figure: %s', ME.message);
end
