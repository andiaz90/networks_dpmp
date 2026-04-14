function obj = smm_obj_for_cmaes(theta, d_mom, W_mat, bas, M_in, opt_in, oo_in)
% SMM_OBJ_FOR_CMAES  Objective function wrapper for cmaes.m
% cmaes requires a named function (not an anonymous handle); this file
% provides that interface and forwards the call to smm_model_moments.
    persistent best_obj_so_far

    % Initialize from previous checkpoint once per MATLAB session so an
    % interrupted run can continue improving from the last saved best value.
    if isempty(best_obj_so_far)
        best_obj_so_far = Inf;
        ckpt_file = fullfile(fileparts(mfilename('fullpath')), 'smm_best_so_far.mat');
        if exist(ckpt_file, 'file')
            try
                tmp = load(ckpt_file, 'smm_best_so_far');
                if isfield(tmp, 'smm_best_so_far') && isfield(tmp.smm_best_so_far, 'obj_best')
                    best_obj_so_far = tmp.smm_best_so_far.obj_best;
                end
            catch
                % Keep Inf if checkpoint loading fails.
            end
        end
    end

    [m_model, info] = smm_model_moments(theta, M_in, opt_in, oo_in, bas);
    if info(1) ~= 0 || any(isnan(m_model))
        obj = 1e8;
        return
    end
    psi = d_mom - m_model;
    obj = psi' * W_mat * psi;

    % Save checkpoint only when a strict improvement is found.
    if obj < best_obj_so_far
        best_obj_so_far = obj;
        theta_best = theta(:);
        nsec = bas.nsec;

        ilabcosts_val  = theta_best(1);
        modepsY        = ones(nsec,1) * theta_best(2);
        modepsM        = ones(nsec,1) * theta_best(3);
        kappaV_val     = exp(theta_best(4));
        rho_om1_val    = theta_best(5);
        sigma_om_val   = theta_best(6);
        rho_tfp1_val   = theta_best(7);
        isigma_tfp_val = theta_best(8:19);
        rho_pvstar_val   = theta_best(20);
        sigma_pvstar_val = theta_best(21);
        rho_xi_val       = theta_best(22);
        sigma_xi_val     = theta_best(23);

        smm_best_so_far = struct( ...
            'theta_best', theta_best, ...
            'obj_best', obj, ...
            'saved_at', datestr(now, 31));

        try
            save(fullfile(fileparts(mfilename('fullpath')), 'smm_best_so_far.mat'), 'smm_best_so_far');
            % Keep this filename updated so main_SOE_gap.m always picks up
            % the best calibration found so far, even mid-optimization.
            save(fullfile(fileparts(mfilename('fullpath')), 'smm_estimates.mat'), ...
                'ilabcosts_val', 'modepsY', 'modepsM', ...
                'kappaV_val', 'rho_om1_val', 'sigma_om_val', 'rho_tfp1_val', ...
                'isigma_tfp_val', 'rho_pvstar_val', 'sigma_pvstar_val', ...
                'rho_xi_val', 'sigma_xi_val');
        catch
            % Optimization should never stop because checkpoint saving failed.
        end
    end
end
