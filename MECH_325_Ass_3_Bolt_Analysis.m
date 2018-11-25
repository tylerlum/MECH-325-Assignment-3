%% MECH 325 Assignment 3
%% Bolt Stress Analysis
%% Note: Currently unset numbers are set to 0

%% Constants From Question
F_t = 37.673; % lbf (force on motor shaft)
t_1 = 0.5;  % in (Thickness of member 1)
t_2 = 1;  % in (Thickness of member 2)

%% Question Parameters %%
d = 1;  % in (Bolt diameter)
L = 3;  % in (Bolt Length)
A_t = 0;  % in^2 (Threaded area)
w = 0.1;  % in (washer thickness)
H = 1;  % in (Nut height)
E = 30 * 10^6;  % psi (Young's Modulus of steel bolt) (Table 8-8)

S_e = 0;  % Endurance Strength
S_ut = 0;  % Ultimate Strength

%% Motor Mount Dimensions
y = 6;  % in (Height of where force is applied)
x = 6;  % in (distance from center to bolt)
m = 6;  % in (length of motor)

%% Bolt calculations
A_d = pi * d^2 / 4;  % in^2 (Bolt shank area) (Table 8-7)
L_t = 2*d + 1/4;  % in (Threaded length) (Eqn 8-13) (L<6 for sure)
l = t_1 + t_2 + 2*w;  % in (grip length) (Table 8-7)
L_d = L - L_t;  % in (Shank/Unthreaded length) (Table 8-7)
l_d = L_d;  % in (Unthreaded length of grip) (Table 8-7)

%% Check if bolt has long enough threaded region
if (L_d >= w + t_1 + t_2)
    error = "L_d is too big, nut can't tighten on member"
end

%% Check if bolt is long enough
if (L <= l + H)
    error = "L too small, can't put on nut"
end

%% Calculate Bolt Tension (diagram in Force Analysis 1.png)
R_b = F_t * (cos(30 * pi / 180) * (x + m / 2) + sin(30 * pi / 180) * y) / (2*x);  % lbf

%% Get stiffness 
k_b = A_d * A_t * E / ((A_t * l_t) + (A_d * l_d));  % lbf/in (bolt stiffness) (Table 8-7)
k_m = 0;  % member stiffness (not sure how to calculate yet)

%% Calculate Safety Factor
F_p = A_t * S_p;
F_i = 0.9 * F_p;


n_fs = S_e * (S_ut - o_i) / ((S_ut * o_a) + S_e * (o_m - o_i));


