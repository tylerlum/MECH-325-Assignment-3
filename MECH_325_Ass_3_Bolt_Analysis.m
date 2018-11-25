%% MECH 325 Assignment 3
%% Bolt Stress Analysis
%% Note: Currently unset numbers are set to 0

%% Constants From Question
F_t = 37.673; % lbf (force on motor shaft)
t_1 = 0.5;  % in (Thickness of member 1)
t_2 = 1;  % in (Thickness of member 2)

%% Question Parameters %% GUESS SAE-5 (0.5 in diameter) (2.5 in length)
d = 0.5;  % in (Bolt diameter)
L = 2.5;  % in (Bolt Length)
A_t = 0.1419;  % in^2 (Threaded area) (Table 8-2) UNC
w = 0.109;  % in (washer thickness) (Table A-32)
H = 7/16;  % in (Nut height) (Table A-31)
E = 30 * 10^6;  % psi (Young's Modulus of steel bolt) (Table 8-8)

S_e = 18.6 * 10^3;  % psi Endurance Strength (Table 8-17)
S_ut = 120 * 10^3;  % psi Ultimate Strength (Table 8-9)
S_p = 85 * 10^3; % psi Proof Strength (Table 8-9)

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
l_t = l - l-d;  % in (threaded length of grip) (Table 8-7)

%% Check if bolt has long enough threaded region
if (L_d > l)
    disp('L_d is too big, nut can not tighten on member')
end

%% Check if bolt is long enough
if (L <= l + H)
    disp('L too small, can not put on nut')
end

%% Calculate Bolt Tension (diagram in Force Analysis 1.png)
P = (cos(30 * pi/180) * (x + m/2) + sin(30 * pi/180) * y) * F_t / (2*x) / 2;  % lbf (divide by 2 at end because two bolts on the high tension side)

%% Get stiffness 
k_b = A_d * A_t * E / ((A_t * l_t) + (A_d * l_d));  % lbf/in (bolt stiffness) (Table 8-7)

%Modulus of Elasticity (psi),Head diameter(in), Bolt diameter(in), Thickness of member (add thickness of washer) (in);
%Assume washer thicnkess is negligable. For K1 and K2 assume head diameter
%is 1.5 times bolt diameter
%See example 8-2
k1 = oneMemberStiffness(E,1.5*d,d,0.5); %Top plate
k2 = oneMemberStiffness(E,1.5*d,d,0.75); %Bottom plate
k3 = oneMemberStiffness(E,(3*d*tan(30 * pi / 180) + d),d,0.25); %Intermediate
k_m = (k1*k2*k3)/(k1*k2 + k1*k3 + k2*k3);  

%% Force calculations
C = k_b / (k_b + k_m);  % fraction of external load carried by bolt (Section 8-7 f)
F_p = A_t * S_p;
F_i = 0.9 * F_p;
% might want to ensure F_m > 0 as per Eqn 8-25
F_m = C*P*k_m/k_b;
if (F_m <= 0)
    disp('F_m is less than 0, separation has occured')
end
%% Stress calculations (Eqns 8-39 to 8-41) We have case P_max = P and P_min = 0, so matches this
o_a = C * P / (2 * A_t);  % psi (alternating stress)
o_i = F_i / A_t;  % psi (preload stress)
o_m = o_a + o_i;  % psi (mean stress)

%% Calculate Safety Factor
n_fs = S_e * (S_ut - o_i) / ((S_ut * o_a) + S_e * (o_m - o_i))  % Goodman safety factor (Eqn 8-38)

%% TODO: Cost calculation (idk how the grade cost equation works)