% Generate a simple example and check solution with Matlab solver
[A, b, c] = RandomLinearProgram(50, 100, 1234);

% My solution
[x, obj] = longStepInteriorPoint(A, b, c, 1000);

% Matlab linear programming solver solution
x_star = linprog(c, [], [], A, b, zeros(100, 1), []);
obj_star = c' * x_star;

% Difference in solution and objectives
fprintf('Error in solution :  %d .\n', norm(x_star - x));
fprintf('Error in objective: %d .\n', norm(obj_star - obj));


