% smm_estimation_cluster.m
%
% Lanzador SMM para el cluster HPC.
% Configura rutas del cluster y llama smm_estimation.m directamente.
%
% Ajustar las tres rutas de la seccion de configuracion segun el servidor.

%% =========================================================================
%%  CONFIGURACION DEL CLUSTER — ajustar segun servidor
%% =========================================================================

% Ruta a Dynare en el cluster
setenv('HPC_DYNARE_PATH', '/repositorio/modules/.local/easybuild/sources/builds/dynare/matlab');

% Ruta a los archivos de datos en el cluster
setenv('HPC_DATA_PATH', '/repositorio/adiaz/networks/data');

% Ruta a utilidades del modelo
addpath('/repositorio/adiaz/networks/utils');

%% =========================================================================
%%  CORRER ESTIMACION
%% =========================================================================
run smm_estimation.m
