% Simplex Examples

A_1 = [1,1,1,0; 2,0.5,0,1];
b_1 = [5; 8];
c_1 = [-4; -2; 0; 0];
% The expected optimal objective is -17.333

[x_1, obj_1, ~, ~, ~] = simplexMethod(A_1, b_1, c_1, 100)

A_2 = [1,1,1,0; 2,0.5,0,1];
b_2 = [10; 12];
c_2 = [-3; -2; 0; 0];

[x_2, obj_2, ~, ~, ~] = simplexMethod(A_2, b_2, c_2, 100)
% The expected optimal objective is -24.667
