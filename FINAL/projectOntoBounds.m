function x_proj = projectOntoBounds(x, l, u)
    % Project vector onto rectangular box [l, u] as defined in (17.52)
    x_proj = min(max(x, l), u);
end