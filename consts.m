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


global A B E
A = [
	-F/V - k_0*exp(-E_R/T), -k_0*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A;
	h*k_0*exp(-E_R/T)/Ro/c_p, -F/V+(k_0*h*V*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A -a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc)))/(V*Ro*c_p)
	];
B = [F_in/V, 0;
	0, -((F_C^b*a*(T - T_Cin)*(b + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro))...
	- (F_C^(b + 1)*a*(T - T_Cin)*((F_C^(b - 1)*a*b)/(2*c_pc*Ro_c) + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro_c))^2)/(V*c_p*Ro)
	];
E = [0, 0;
	F_in/V, a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
	];