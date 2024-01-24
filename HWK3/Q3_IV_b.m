% Define initial guess, tolerance, and maximum number of iterations
x0_2 = [3; 4];
tol = 1e-6;  % Tolerance for convergence
max_iter = 100; % Maximum number of iterations

% Apply the VN method to Fenton's function
x_star_VN_2 = VN(@fentonfgH, x0_2, tol, max_iter);
format shortE

fprintf('VN initialized at [3, 4] ended at x = [%e, %e]\n', x_star_VN_2(1), x_star_VN_2(2));

% Apply the NMHM method to Fenton's function
x_star_NMHM_2 = NMHM(@fentonfgH, x0_2, tol, max_iter);
fprintf('NMHM initialized at [3, 4] found a solution at x = [%f, %f]\n', x_star_NMHM_2(1), x_star_NMHM_2(2));
