% Temporary wrapper to run analysis
cd('c:\github\networks_dpmp_dme\DPMP-DME-Modelos\Model\modelo_chile');
load('model_output_IOSOE_ex2.mat');
load('params_val.mat');
analyze_services_demand;
