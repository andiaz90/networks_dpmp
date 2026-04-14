%% Test just the FlexNoIO config to debug the issue

clear all; close all; clc;

restoredefaultpath;
set(0,'DefaultLineLineWidth',2);

this_dir = cd;
addpath(fullfile(this_dir, '..', 'utils'));
addpath('C:\Program Files\Dynare\6.4\matlab');

nsec = 12;

% Load baseline parameters
load('baseline_params.mat', 'modbeta_full', 'modbeta_diag', 'kappa_flex', 'kappa_sticky', 'alpha', 'modalphaV', 'modvarrho', 'omegaX_val', 'names');
load('IO_data.mat', 'spend_good', 'spend_serv');

fprintf('\n=== Testing FlexNoIO Configuration ===\n\n');

%% Set up FlexNoIO config
config_name = 'FlexNoIO';
config_desc = 'Flexible Price + Diagonal I-O + No Labor Costs';

fprintf('Config: %s\n', config_desc);
fprintf('Beta matrix: diagonal I-O\n');
fprintf('Kappa: flexible prices\n');
fprintf('Labor costs: DISABLED\n\n');

% Assign parameters to workspace for steady state computation
modbeta  = modbeta_diag;
modkappa = kappa_flex;
modalpha = alpha;
modcl    = zeros(nsec,1);  % NO labor costs
modclneg = zeros(nsec,1);
modcm    = zeros(nsec,1);
modalphaV_cfg = modalphaV;
modvarrho_cfg = modvarrho;
omegaX_cfg = omegaX_val;

assignin('base', 'nsec',      nsec);
assignin('base', 'modbeta',   modbeta);
assignin('base', 'modkappa',  modkappa);
assignin('base', 'modalpha',  modalpha);
assignin('base', 'modalphaV', modalphaV_cfg);
assignin('base', 'modcl',     modcl);
assignin('base', 'modclneg',  modclneg);
assignin('base', 'modcm',     modcm);
assignin('base', 'modvarrho', modvarrho_cfg);
assignin('base', 'omegaX',    omegaX_cfg);

% Load scalar parameters
load('scalar_params.mat');

fprintf('Attempting steady state computation for FlexNoIO...\n\n');

% Check if modbeta is diagonal
fprintf('Beta matrix properties:\n');
fprintf('  Shape: %d x %d\n', size(modbeta, 1), size(modbeta, 2));
fprintf('  Trace: %.6f\n', trace(modbeta));
fprintf('  Off-diagonal max: %.6e\n', max(max(abs(modbeta - diag(diag(modbeta))))));
fprintf('  Rank: %d\n', rank(modbeta));
fprintf('\n');

% Set up kappa
fprintf('Kappa properties:\n');
fprintf('  Mean: %.6f\n', mean(modkappa));
fprintf('  Min: %.6f\n', min(modkappa));
fprintf('  Max: %.6f\n', max(modkappa));
fprintf('  Flexible prices? %d\n', max(modkappa) < 1e-10);
fprintf('\n');

% Now try to save and run Dynare
try
    % Get steady state values needed for Dynare initialization
    fprintf('Computing steady state...\n');
    diary('test_flexnoio_log.txt');
    
    % Call Dynare
    fprintf('Running Dynare for FlexNoIO config...\n');
    dynare NK_IOSOE_tot.mod noclearall;
    
    diary off;
    fprintf('\n✓ SUCCESS: FlexNoIO config completed\n');
    
catch ME
    diary off;
    fprintf('\n✗ ERROR in FlexNoIO config:\n');
    fprintf('Message: %s\n', ME.message);
    fprintf('\nStack trace:\n%s\n', ME.getReport('extended'));
    fprintf('\nDebugging info:\n');
    fprintf('  Last error: %s\n', lasterr);
    
    % Check workspace for relevant variables
    fprintf('\n  Workspace check:\n');
    fprintf('  - modkappa exists: %d\n', exist('modkappa', 'var') > 0);
    if exist('modkappa', 'var')
        fprintf('    modkappa flexible? %d\n', max(modkappa) < 1e-10);
    end
    fprintf('  - modbeta exists: %d\n', exist('modbeta', 'var') > 0);
    if exist('modbeta', 'var')
        fprintf('    modbeta rank: %d (expected %d)\n', rank(modbeta), nsec);
    end
end
