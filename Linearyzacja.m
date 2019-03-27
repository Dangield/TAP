close all
clear variables
global Ro Ro_c c_p c_pc k_0 E_R h a b;
global V F_in F C_Ain F_C T_in T_Cin;

%Sta³e
Ro = 1e6;
Ro_c = 1e6;
c_p = 1;
c_pc = 1;
k_0 = 1e10;
E_R = 8330.1;
h = 130e6;
a = 1.678e6;
b = 0.5;

% Wartoœci steruj¹ce
V = 1;
F_in = 1;
F = 1;
% Sterowanie
C_Ain = 2;
F_C = 15;
% Zak³ócenie
T_in = 323;
T_Cin = 365;
% Wyjœcia
C_A = 0.26;
T = 393.9;

% linearyzacja:
global A B E
A = [
	-F/V - k_0*exp(-E_R/T), -k_0*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A;
	h*k_0*exp(-E_R/T)/Ro/c_p, -F/V-a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
	]/60;
eig(A)
B = [F_in/V, 0;
	0, -((F_C^b*a*(T - T_Cin)*(b + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro))...
	- (F_C^(b + 1)*a*(T - T_Cin)*((F_C^(b - 1)*a*b)/(2*c_pc*Ro_c) + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro_c))^2)/(V*c_p*Ro)
	]/60;
E = [0, 0;
	F_in/V, a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
	]/60;

% Zmienne do symulacji
t0 = 0;
tfinal = 500;
x0 = [C_A; T];
coef_vals = [0.6, 0.9, 0.95, 1, 1.05, 1.1, 1.4];

% Zmienne do plotowania
labels = [];
color_formats = [[1 0 0]; [0 1 0]; [0 0 1]; [0 0.7 0.7]; [0.7 0 0.7]; [0.7 0.7 0]; [0 0 0]];

figure(1)
for i=1:7
	
	c0 = [C_Ain*coef_vals(i); F_C];
	d0 = [T_in; T_Cin];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	plot(t, x(:, 1), 'Color', color_formats(i, :))
	hold on
	
	labels = [labels, "Model nieliniowy, C_Ain zmnienione " + coef_vals(i) + " razy"];
end

lin_x0 = x0 - [C_A; T];

for i=1:7
	c0 = [C_Ain*coef_vals(i); F_C];
	d0 = [T_in; T_Cin];
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 1) = lin_x(:, 1) + C_A;
	
	plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	hold on
	
	labels = [labels, "Model liniowy, C_Ain zmnienione " + coef_vals(i) + " razy"];
end
legend(labels)
title("C_A")

figure
for i=1:7
	c0 = [C_Ain*coef_vals(i); F_C];
	d0 = [T_in; T_Cin];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	plot(t, x(:, 2), 'Color',  color_formats(i, :))
	hold on
end

for i=1:7
	coef = coef_vals(i);
	c0 = [C_Ain*coef; F_C];
	d0 = [T_in; T_Cin];
	
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 2) = lin_x(:, 2) + T;
	
	plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
	hold on
end
title("T")
legend(labels)


function dx = reactor(t, x, input, disturbance)
	global Ro Ro_c c_p c_pc k_0 E_R h a b;
	global V F_in F;

	d_Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))/60;
	d_T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))/60;

	dx = [d_Ca; d_T];
end

function dx = lin_reactor(t, x, input, disturbance)
	global A B E;

	dx = A*x + B*input + E*disturbance;
end