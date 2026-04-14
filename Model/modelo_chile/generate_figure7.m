% generate_figure7.m
% Complete workflow to generate Figure 7: Manufacturing TFP Shock Effects
% 
% This script:
% 1. Runs the model with manufacturing TFP shock (exercise 2)
% 2. Generates Figure 7 showing sectoral and macro effects
% 3. Saves high-resolution figures for the paper

fprintf('\n');
fprintf('============================================================\n');
fprintf('   FIGURE 7 GENERATION: MANUFACTURING TFP SHOCK ANALYSIS   \n');
fprintf('============================================================\n');
fprintf('\n');

%% Step 1: Check if model results already exist
fprintf('STEP 1: Checking for existing model results...\n');

if exist('model_output_IOSOE_ex2.mat', 'file')
    fprintf('  → Found existing model output: model_output_IOSOE_ex2.mat\n');
    fprintf('  → Loading results...\n');
    
    load('model_output_IOSOE_ex2.mat');
    
    if exist('oo_', 'var') && isfield(oo_, 'irfs')
        fprintf('  ✓ Model results loaded successfully\n');
        fprintf('  ✓ IRFs available: %d variables\n', length(fieldnames(oo_.irfs)));
        use_existing = true;
    else
        fprintf('  ✗ Loaded file incomplete, will re-run model\n');
        use_existing = false;
    end
else
    fprintf('  → Model output not found, will run model\n');
    use_existing = false;
end

%% Step 2: Run model if needed
if ~use_existing
    fprintf('\nSTEP 2: Running manufacturing TFP shock simulation...\n');
    fprintf('  This may take 5-10 minutes...\n\n');
    
    % Run exercise 2
    run('main_SOE_exercise2.m');
    
    fprintf('\n  ✓ Model simulation complete\n');
else
    fprintf('\nSTEP 2: Skipping model run (using existing results)\n');
end

%% Step 3: Generate Figure 7
fprintf('\nSTEP 3: Generating Figure 7...\n');

% Run the plotting script
run('plot_figure7_manufacturing_shock.m');

fprintf('\n✓ All figures generated successfully!\n');

%% Step 4: Display file locations
fprintf('\n============================================================\n');
fprintf('   FIGURE 7 - FILES GENERATED                              \n');
fprintf('============================================================\n');
fprintf('\n');
fprintf('Figure 7a (Sectoral Effects):\n');
fprintf('  → figure7a_sectoral_effects.png\n');
fprintf('\n');
fprintf('Figure 7b (Macro Effects):\n');
fprintf('  → figure7b_macro_effects.png\n');
fprintf('\n');
fprintf('Figure 7 (Combined View):\n');
fprintf('  → figure7_manufacturing_shock_complete.png\n');
fprintf('\n');
fprintf('Model Output:\n');
fprintf('  → model_output_IOSOE_ex2.mat\n');
fprintf('\n');

%% Step 5: Summary instructions
fprintf('============================================================\n');
fprintf('   NEXT STEPS                                              \n');
fprintf('============================================================\n');
fprintf('\n');
fprintf('To insert figures in your LaTeX paper:\n');
fprintf('\n');
fprintf('\\begin{figure}[H]\n');
fprintf('  \\centering\n');
fprintf('  \\includegraphics[width=\\textwidth]{figure7_manufacturing_shock_complete.png}\n');
fprintf('  \\caption{Impact of Manufacturing TFP Shock on Sectoral and Macroeconomic Variables}\n');
fprintf('  \\label{fig:manufacturing_shock}\n');
fprintf('\\end{figure}\n');
fprintf('\n');
fprintf('Alternatively, use separate panels:\n');
fprintf('- figure7a_sectoral_effects.png for sectoral variables\n');
fprintf('- figure7b_macro_effects.png for macro/OE variables\n');
fprintf('\n');
fprintf('============================================================\n');
fprintf('\n');
