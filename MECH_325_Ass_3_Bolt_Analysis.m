%% MECH 325 Assignment 3
%% Bolt Stress Analysis

%% Constants From Question
F_t = 0; % lbf (force on motor shaft)
t_1 = 0.5;  % in (Thickness of member 1)
t_2 = 1;  % in (Thickness of member 2)

%% Bolt sizing
d = 1;  % in (Bolt diameter)
L = 3;  % in (Bolt Length)
A_d = pi * d^2 / 4;  % in^2 (Bolt shank area)
A_t = 0;  % in^2 (Threaded area)
L_t = 2*d + 1/4;  % in (8-13) (Threaded length)

%% Washer and nut
w = 0.1;  % in (Washer thickness)
H = 1;  % in (Nut height)

l = t_1 + t_2 + 2*w;
L_d = L - L_t;

%% Motor Mount Dimensions
y = 6;
m = 6;
x = 6;

%% Check if bolt has long enough threaded region
if (L_d >= w + t_1 + t_2)
    error = "L_d is too small, nut can't tighten on member"
end

%% Check if bolt is long enough
if (L <= l = H)
    error = "L too small, can't put on nut"
end

    
    



%% Calculate Bolt Tension 
R_b = F_t * (cos(30 * pi / 180) * (x + m / 2) + sin(30 * pi / 180))/ (2 * x);

%% Get stiffness 
k_b = A_d * A_t * E / ((A_t * l_t) + (A_d * l_d));
k_m = 0;

%% Calculate Safety Factor
F_p = A_t * S_p;
F_i = 0.9 * F_p;


n_fs = S_e * (S_ut - o_i) / ((S_ut * o_a) + S_e * (o_m - o_i));


