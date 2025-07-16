function [U_mean, norm_V]= karcher_mean_grassmann(projector_matrices, epsilon, tol, max_iter,k)
    % Karcher Mean Computation on Grassmannian using Projector Matrices
    %
    % Inputs:
    %   projector_matrices - Cell array of k (NxN) projector matrices (q_i)
    %   epsilon - Step size for the update (typically 0.5)
    %   tol - Tolerance for convergence (e.g., 1e-6)
    %   max_iter - Maximum number of iterations
    %
    % Output:
    %   U_mean - The Karcher mean projector matrix (u_j)

    % Initialization
     % k Number of points
    U_mean = projector_matrices{1}; % Initial mean (randomly chosen point)
    nq = size(U_mean, 1);
    
    % Identity matrix
    
    I = [eye(nq) zeros(nq) zeros(nq) zeros(nq)];

    % Iterative procedure
    for iter = 1:max_iter
        % Compute tangent vectors v_i using the inverse exponential map
        V = [];%zeros(N, N, k); % Tangent vectors v_i
        for i = 1:k
            Q = projector_matrices{i};
diff_matrix = timesQ((I - 2 * Q), (I - 2 * U_mean));
     log_mat = 0.5*logm(RealQuaternion(diff_matrix));
             log_mat=QuaternionReal(log_mat);
           % log_mat = 0.5 * logm(I - 2 * Q) - 0.5 * logm(I - 2 * U_mean);
            V(:, :, i) = timesQ(U_mean, log_mat) - timesQ(log_mat , U_mean);
        end
        
        % Average tangent vector
        V_avg = [zeros(nq), zeros(nq), zeros(nq), zeros(nq)];
        for i = 1:k
            V_avg = V_avg + V(:, :, i);
        end
        V_avg = V_avg / k;
        
        % Check convergence
        norm_V = norm(V_avg, 'fro');
        if norm_V < tol
            fprintf('Converged at iteration %d with norm %.6f\n', iter, norm_V);
            break;
        end
        
        % Update using the exponential map
        exp_term = expm(RealQuaternion(epsilon * timesQ(V_avg, U_mean) - timesQ(U_mean, epsilon * V_avg)));

        exp_term=QuaternionReal(exp_term);

        exp_termn = expm(RealQuaternion(-epsilon * timesQ(V_avg, U_mean) + timesQ(U_mean, epsilon * V_avg)));

        exp_termn=QuaternionReal(exp_termn);

        U_mean = timesQ(timesQ(exp_term, U_mean), exp_termn);
        % Symmetrize to ensure U_mean remains a projector matrix
        U_mean = (U_mean + transQ(U_mean)) / 2;
    end
    
    if iter == max_iter
        fprintf('Reached maximum iterations without full convergence.\n');
    end
end