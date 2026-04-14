%% submit_smm_cluster.m
%
% Envia el trabajo de estimacion SMM al cluster HPC.
%
% FLUJO:
%   1. Verifica prerequisitos en la maquina local (compila el modelo si hace
%      falta, genera data_moments_chile.mat si no existe).
%   2. Guarda el workspace de Dynare (M_, oo_, options_) en
%      dynare_workspace.mat.
%   3. Envia smm_estimation_cluster.m al cluster con todos los archivos
%      de datos y funciones necesarios adjuntos.
%
% PREREQUISITOS LOCALES:
%   - MATLAB + Dynare instalados localmente.
%   - El modelo NK_IOSOE_lev.mod debe estar compilado al menos una vez
%     (la primera ejecucion de main_SOE.m lo hace). Si M_/oo_/options_
%     no estan en el workspace, este script corre main_SOE.m.
%   - compute_data_moments.m debe haberse corrido (genera
%     data_moments_chile.mat).
%
% INSTRUCCIONES DE USO:
%   1. Ajustar REPO_FOLDER con la ruta en el repositorio del cluster.
%   2. Correr este script desde la carpeta Model/modelo_chile/.
%   3. Tras enviar el job, monitorear con:
%          >> myjob.State
%      Para recuperar resultados:
%          >> load(myjob)
%          >> diary(myjob)
% =========================================================================

%% =========================================================================
%%  CONFIGURACION DEL USUARIO — ajustar antes de enviar
%% =========================================================================

% Perfil del cluster configurado en MATLAB Parallel Computing Toolbox
CLUSTER_PROFILE = 'HPC R2024b';

% Numero de workers para el pool paralelo.
CORES = 40;

% Carpeta de trabajo en el repositorio del cluster.
% Debe existir previamente. Modificar segun usuario y proyecto.
REPO_FOLDER = '/repositorio/UsuarioBanco/DPMP-DME/';

%% =========================================================================
%%  1. VERIFICAR PREREQUISITOS EN LA MAQUINA LOCAL
%% =========================================================================
fprintf('\n=============================================================\n');
fprintf('  PREPARACION DEL JOB SMM PARA EL CLUSTER\n');
fprintf('=============================================================\n\n');

% El modelo se compilara directamente en el cluster con main_SOE.m.
% Solo se necesita data_moments_chile.mat generado localmente.
if ~exist('data_moments_chile.mat','file')
    fprintf('[1] data_moments_chile.mat no encontrado.\n');
    fprintf('    Corriendo compute_data_moments.m...\n');
    try
        run compute_data_moments.m
    catch ME_dm
        error('compute_data_moments.m fallo: %s', ME_dm.message);
    end
    if ~exist('data_moments_chile.mat','file')
        error('compute_data_moments.m no genero data_moments_chile.mat.');
    end
    fprintf('[1] data_moments_chile.mat generado OK.\n');
else
    fprintf('[1] data_moments_chile.mat presente.\n');
end

%% =========================================================================
%%  2. ENVIAR JOB AL CLUSTER
%% =========================================================================

% Solo se adjunta data_moments_chile.mat (archivo pequeno generado localmente).
% El resto de los archivos (main_SOE.m, NK_IOSOE_lev.mod, smm_model_moments.m,
% steady_ntwsoe*.m, etc.) deben estar en REPO_FOLDER en el servidor.
attached_files = {
    'data_moments_chile.mat'     ... % momentos de datos de Chile
};

% Verificar existencia de cada archivo adjunto
for k = 1:numel(attached_files)
    if ~exist(attached_files{k}, 'file')
        warning('Archivo adjunto no encontrado: %s', attached_files{k});
    end
end

% Inicializar conexion con el cluster
fprintf('[2] Conectando con el cluster "%s"...\n', CLUSTER_PROFILE);
clust = parcluster(CLUSTER_PROFILE);

fprintf('    Enviando smm_estimation_cluster.m al cluster...\n');
fprintf('    Workers:       %d\n', CORES);
fprintf('    Carpeta repo:  %s\n\n', REPO_FOLDER);

myjob = batch(clust, 'smm_estimation_cluster', ...
    'Pool',              CORES, ...
    'CurrentFolder',     REPO_FOLDER, ...
    'AutoAddClientPath', false, ...
    'AttachedFiles',     attached_files);

%% =========================================================================
%%  4. INFORMACION POST-ENVIO
%% =========================================================================
fprintf('=============================================================\n');
fprintf('  JOB ENVIADO AL CLUSTER\n');
fprintf('=============================================================\n\n');
fprintf('  Job ID:    %s\n',   myjob.ID);
fprintf('  Estado:    %s\n\n', myjob.State);
fprintf('Comandos utiles:\n');
fprintf('  Estado del job:          myjob.State\n');
fprintf('  Esperar a que termine:   wait(myjob)\n');
fprintf('  Cargar resultados:       load(myjob)\n');
fprintf('  Ver salida de consola:   diary(myjob)\n');
fprintf('  Cancelar el job:         cancel(myjob)\n');
fprintf('  Borrar el job:           delete(myjob)\n\n');
fprintf('Los archivos de resultados se guardan en:\n');
fprintf('  %ssmm_results_cluster.mat\n',   REPO_FOLDER);
fprintf('  %ssmm_estimates_cluster.mat\n', REPO_FOLDER);
