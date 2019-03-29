global Ro Ro_c c_p c_pc k_0 E_R h a b;
global V F_in F C_Ain F_C T_in T_Cin C_A T;

%% Stale
Ro = 1e6;
Ro_c = 1e6;
c_p = 1;
c_pc = 1;
k_0 = 1e10;
E_R = 8330.1;
h = 130e6;
a = 1.678e6;
b = 0.5;

%% Wartosci sterujace
V = 1;
F_in = 1;
F = 1;

% Sterowanie
C_Ain = 2;
F_C = 15;

% Zaklocenie
T_in = 323;
T_Cin = 365;

% Wyjscia
C_A = 0.2646;
T = 393.9531;