function [f, g, H] = fentonfgH(x)
    x1 = x(1);
    x2 = x(2);

    % Function value
    f = (12 + x1^2 + (1 + x2^2) / x1^2 + 1 / (x1^2 * x2^2) + 100 / (x1 * x2)^4) / 10;

    % Gradient
    dfdx1 = (2*x1 - 2*(1 + x2^2)/x1^3 - 2/(x1^3 * x2^2) - 400/(x1^5 * x2^4)) / 10;
    dfdx2 = (2*x2/x1^2 - 2/(x1^2 * x2^3) - 400/(x1^4 * x2^5)) / 10;
    g = [dfdx1; dfdx2];

    % Hessian
    d2fdx1x1 = (2 + 6*(1 + x2^2)/x1^4 + 6/(x1^4 * x2^2) + 2000/(x1^6 * x2^4)) / 10;
    d2fdx1x2 = (-4*x2/x1^3 + 4/(x1^3 * x2^3) + 1600/(x1^5 * x2^5)) / 10;
    d2fdx2x2 = (2/x1^2 + 6/(x1^2 * x2^4) + 2000/(x1^4 * x2^6)) / 10;
    H = [d2fdx1x1, d2fdx1x2; d2fdx1x2, d2fdx2x2];
end