% Define initial guess, tolerance, and maximum number of iterations
x0_1 = [3; 2];
tol = 1e-6;  % Tolerance for convergence
max_iter = 100; % Maximum number of iterations

% Apply the VN method to Fenton's function
x_star_VN_1 = VN(@fentonfgH, x0_1, tol, max_iter);
fprintf('VN initialized at [3, 2] found a solution at x = [%f, %f]\n', x_star_VN_1(1), x_star_VN_1(2));

% Apply the NMHM method to Fenton's function
x_star_NMHM_1 = NMHM(@fentonfgH, x0_1, tol, max_iter);
fprintf('NMHM initialized at [3, 2] found a solution at x = [%f, %f]\n', x_star_NMHM_1(1), x_star_NMHM_1(2));

