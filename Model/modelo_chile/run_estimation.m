clear
close all
clc
restoredefaultpath
set(0, 'DefaultLineLineWidth', 2);
warning off
set(0,'DefaultFigureWindowStyle','docked')
restoredefaultpath

% Change this path accordingly
restoredefaultpath
addpath N:\reallocation\FGI\4.5.6\matlab
addpath N:\FGI_Replication\Model\utils



delete lsq*

nsec_val = 14;

filename_industries='Stata_to_excel_few_industries.xls';
filename_io = 'IOmat_few_industries.csv';
filename_fpa='fpa_vector_few_industries.csv';

isec_val = 1:nsec_val;
idropsec = 0;



%% READ DATA
varNames = {'kk','BEA_line','BEACode','BEAName','Nameshort',...
    'da','dp','dy','dl','dri'...
    'share','industrytype','spend_good','spend_serv','alpha'};
varTypes = {'double','double','char','char', 'char', ...
    'double','double','double','double','double',...
    'double','categorical','double','double','double'};
opts = spreadsheetImportOptions('NumVariables',15,...
    'VariableNames',varNames,'VariableTypes',varTypes,'DataRange','A2');
alldata = readtable(filename_industries,opts);
names = table2array(alldata(:,5));
share = table2array(alldata(:,11));
tfpshock = table2array(alldata(:,6))/100;
p_d = table2array(alldata(:,7));
y_d = table2array(alldata(:,8));
l_d = table2array(alldata(:,9));
ri_d = table2array(alldata(:,10));
spend_good = table2array(alldata(:,13));
spend_serv = table2array(alldata(:,14));
alpha = table2array(alldata(:,15));



%% DEFINITIONS

fid = fopen('definition_block_nsec.mod', 'w');
fprintf(fid, '%s\r\n', ['@#define nsec = '  num2str(nsec_val)] ) ;
type definition_block_nsec.mod;
fclose(fid);

fid = fopen('definition_block_io.mod', 'w');
fprintf(fid, '%s\r\n', ['@#define io_val = 1 ' ]) ;
type definition_block_io.mod;
fclose(fid);

fid = fopen('definition_block_lab.mod', 'w');
fprintf(fid, '%s\r\n', ['@#define I_lab_sec = 1 ' ]) ;
fprintf(fid, '%s\r\n', ['@#define I_hiring_cost = 1 ' ]) ;
fprintf(fid, '%s\r\n', ['@#define I_asym_cost = 0 ' ]) ;
fprintf(fid, '%s\r\n', ['@#define I_asym_cost2 = 0 ' ]) ;
fprintf(fid, '%s\r\n', ['@#define I_sym_cost = 0 ' ]) ;
type definition_block_lab.mod;
fclose(fid);


%% SOLUTION BLOCK

fid = fopen('solution_block.mod', 'w');

fprintf(fid, '%s\r\n', ['shocks; ' ]) ;
fprintf(fid, '%s\r\n', ['var eps_om; ' ]) ;
fprintf(fid, '%s\r\n', ['periods 1;']);
fprintf(fid, '%s\r\n', ['values (1);']);

fprintf(fid, '%s\r\n', ['var eps_i; ' ]) ;
fprintf(fid, '%s\r\n', ['periods 1;']);
fprintf(fid, '%s\r\n', ['values (1);']);

fprintf(fid, '%s\r\n', ['var epschi; ' ]) ;
fprintf(fid, '%s\r\n', ['periods 1;']);
fprintf(fid, '%s\r\n', ['values (1);']);

fprintf(fid, '%s\r\n', ['@#for z in 1:nsec' ]) ;
fprintf(fid, '%s\r\n', ['   var epsA_@{z};' ]) ;
fprintf(fid, '%s\r\n', ['   periods 1;' ]) ;
fprintf(fid, '%s\r\n', ['   values (tfpshock(@{z}));' ]) ;
fprintf(fid, '%s\r\n', ['@#endfor' ]) ;
fprintf(fid, '%s\r\n', ['end;']);


fprintf(fid, '%s\r\n', ['perfect_foresight_setup(periods=150);']);
fprintf(fid, '%s\r\n', ['perfect_foresight_solver;']);
type solution_block.mod;
fclose(fid);

%% SET PARAMETERS

betaio = readmatrix(filename_io);
betaio = betaio(isec_val,isec_val);
betax = betaio./sum(betaio,1);
betaxt = betax';
modbeta = betaxt;

modalpha = alpha ; 

beta_val = 0.995;

fpa = readmatrix(filename_fpa);
fpa = fpa(isec_val);
theta = (1-fpa).^3; % convert to quarterly prob(cannot adjust)
kappa = theta*(10-1)./((1-theta).*(1-theta*beta_val)); % convert to rotemberg
modkappa = kappa;

goods = zeros(nsec_val,1);
goods(spend_good>spend_serv) = 1;
services = zeros(nsec_val,1);
services(spend_serv>spend_good) = 1;
other = zeros(nsec_val,1);
other(goods+services==0) = 1;
goods = logical(goods);
services = logical(services);
other = logical(other);

gammag = spend_good;
gammag = gammag/sum(gammag);
modgammag = gammag;

gammas = spend_serv;
gammas = gammas/sum(gammas);
modgammas = gammas;

% Calibrated parameters
sigma_i_val = 0;
phi_val = 1.5;
rhoi_val=0.0;

rho_om1_val = 0.975;
rho_om2_val = 0.00;

sigma_tfp = 1;
rho_tfp1_val = 0.95;
rho_tfp2_val = 0.00;
[yyar2, maxyyrel, tt]=sim_ar2([rho_tfp1_val rho_tfp2_val 1],20,0);
tfpshock = tfpshock/mean(yyar2(1:8)); %rescaling

rho_val = 0.95;
sigma_L_agg_val = 0.04;

modcl = zeros(nsec_val,1);
modcl(:) = 50;

modclneg = zeros(nsec_val,1);
modclneg(:) = 0;
modcm = zeros(nsec_val,1);

modepsM = ones(nsec_val,1);
modepsM(:) = 0.1;

modepsY = ones(nsec_val,1);
modepsY(:) = 0.8;

isigma_tfp_val=1;

rhoi_val=0.0;

ombar_val = 0.31;
sigma_om_val = 0.045;

ilabcosts_val = 1;
rhoirule_val = 0.0;

%% CALIBRATION V. ESTIMATION
isigma_tfp_val=1;
sigma_L_agg_val = 0.0;

modcl = ones(nsec_val,1)*20;
modclneg = ones(nsec_val,1)*0;
modepsM = ones(nsec_val,1)*0.0906;
modepsY = ones(nsec_val,1)*0.6081;
sigma_L_agg_val = 0.0509;

save('params_val',...
    'modgammag','modgammas','modalpha','modkappa','goods','services',...
    'modbeta','modepsM','modepsY','tfpshock',...
    'rho_val','rhoi_val','rho_om1_val','rho_om2_val','rho_tfp1_val','rho_tfp2_val',...
    'modcl','modclneg','modcm','rhoirule_val','phi_val','ombar_val',...
    'sigma_L_agg_val','sigma_i_val','sigma_om_val','isigma_tfp_val',...
    'beta_val','ilabcosts_val');

dynare NK_IO_8Mar23 noclearall

numit = 1 ;
save numit numit


%% Call function
options.FunctionTolerance = 1e-5;
options.MaxIterations = 100;

param_name = {'modcl','modepsM','modepsY','sigma_L_agg_val'};
[xxx,RESNORM,RESIDUAL,EXITFLAG,OUTPUT,LAMBDA,JACOBIAN] = ...
lsqnonlin(@(x)loss_model(x,param_name,oo_,M_,options_,nsec_val,names,p_d,y_d,l_d,goods,services,share),...
    [ 18.8    0.13    0.82  0.09],...
    [ 1        0.01  0.01  0.00 ],...
    [ 1000     0.99  0.99  0.30 ],...
    options)  ;

save estimation_results_jacobian

load estimation_results_jacobian

VVV = inv(JACOBIAN'*JACOBIAN);

SE = diag(VVV).^0.5;

format short
disp('Coefficients')
fprintf(' %2.2f \n', full(xxx)); 
disp('---')
disp('S.E.')
fprintf(' %2.2f \n', full(SE)); 


