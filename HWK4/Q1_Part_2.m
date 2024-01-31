% Initializations
n_values = 200:100:20000;
a = zeros(size(n_values));
b = zeros(size(n_values));

% Measure total time of experiments when start
start = tic;

for i = 1:length(n_values)
    n = n_values(i);
    fprintf('n = %d\n', n);

    x = 2 * ones(n, 1); % Point of interest
    H = rosenbrocknHessian(x);

    % Timing for sparse format
    H_sparse = sparse(H);
    tic;
    x_sparse = H_sparse \ ones(n, 1);
    a(i) = toc;

    % Timing for dense format
    tic;
    x_dense = H \ ones(n, 1);
    b(i) = toc;

    % Stop iterations when total time elapsed is over 3 minutes
    if toc(start) > 180
        fprintf('Total script time exceeded 3 minutes, breaking loop at n = %d\n', n);
        break;
    end
end

% Plot the results
figure;
ns = n_values(1:i);
as = a(1:i);
bs = b(1:i);

% Sparse
% check approximate power law
subplot(1, 2, 1);
p = polyfit(log(ns),log(as),1); % y = C x^k
C_1 = p(1);
k_1 = exp(p(2));
plot(ns, as, 'b-', ...
    ns, k_1*ns.^C_1, 'k-.','LineWidth', 1.5);
% plotting against a straight line
legend('Sparse a(n)', sprintf('y = k * n^{%.3g}', C_1));
xlabel('Order n of Rosenbrock function');
ylabel('Compute time (s)');
title('Solving a sparse Hessian');

% Dense
% check approximate power law
subplot(1, 2, 2);
p = polyfit(log(ns),log(bs),1); % y = C x^k
C_2 = p(1);
k_2 = exp(p(2));
plot(ns, bs, 'r-', ...
    ns, k_2*ns.^C_2, 'k-.','LineWidth', 1.5);
% plotting against a cubic
legend('Dense b(n)', sprintf('y = k * n^{%.3g}', C_2));
xlabel('Order n of Rosenbrock function');
ylabel('Compute time (s)');
title('Solving a dense Hessian');
