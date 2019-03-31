close all

% Zmienne do plotowania
labels = [];
color_formats = [[1 0 0]; [0 1 0]; [0 0 1]; [0 0.7 0.7]; [0.7 0 0.7]; [0.7 0.7 0]; [0 0 0]];

% Zmienne do symulacji
t0 = 0;
tfinal = 20;
x0 = [C_A; T];
coef_vals = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6];

%% zmiana 1. sterowania

figure(1)
figure(2)
for i=1:7
	
	c0 = [C_Ain*coef_vals(i); F_C];
	d0 = [T_in; T_Cin];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	figure(1)
	hold on
	plot(t, x(:, 1), 'Color', color_formats(i, :))
	figure(2)
	hold on
	plot(t, x(:, 2), 'Color', color_formats(i, :))
	
	labels = [labels, "Model nieliniowy, C_{Ain} zmnienione " + coef_vals(i) + " razy"];
end

lin_x0 = x0 - [C_A; T];

for i=1:7
	c0 = [C_Ain*coef_vals(i); F_C];
	d0 = [T_in; T_Cin];
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 1) = lin_x(:, 1) + C_A;
	lin_x(:, 2) = lin_x(:, 2) + T;
	
	figure(1)
	hold on
	plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(2)
	hold on
	plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
	
	labels = [labels, "Model liniowy, C_{Ain} zmnienione " + coef_vals(i) + " razy"];
end
figure(1)
% legend(labels)
title("C_A(C_{Ain})")

figure(2)
title("T(C_{Ain})")

%% zmiana 2. sterowania
labels = [];

figure(3)
figure(4)
for i=1:7
	
	c0 = [C_Ain; F_C*coef_vals(i)];
	d0 = [T_in; T_Cin];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	figure(3)
	hold on
	plot(t, x(:, 1), 'Color', color_formats(i, :))
	figure(4)
	hold on
	plot(t, x(:, 2), 'Color', color_formats(i, :))
	
	labels = [labels, "Model nieliniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end

lin_x0 = x0 - [C_A; T];

for i=1:7
	c0 = [C_Ain; F_C*coef_vals(i)];
	d0 = [T_in; T_Cin];
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 1) = lin_x(:, 1) + C_A;
	lin_x(:, 2) = lin_x(:, 2) + T;
	
	figure(3)
	hold on
	plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(4)
	hold on
	plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
	
	labels = [labels, "Model liniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end
figure(3)
% legend(labels)
title("C_A(F_C)")

figure(4)
title("T(F_C)")

%% Zmiana 1. zakłócenia:
labels = [];

figure(5)
figure(6)
for i=1:7
	
	c0 = [C_Ain; F_C];
	d0 = [T_in*coef_vals(i); T_Cin];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	figure(5)
	hold on
	plot(t, x(:, 1), 'Color', color_formats(i, :))
	figure(6)
	hold on
	plot(t, x(:, 2), 'Color', color_formats(i, :))
	
	labels = [labels, "Model nieliniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end

lin_x0 = x0 - [C_A; T];

for i=1:7
	c0 = [C_Ain; F_C];
	d0 = [T_in*coef_vals(i); T_Cin];
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 1) = lin_x(:, 1) + C_A;
	lin_x(:, 2) = lin_x(:, 2) + T;
	
	figure(5)
	hold on
	plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(6)
	hold on
	plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
	
	labels = [labels, "Model liniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end
figure(5)
% legend(labels)
title("C_A(T_{in})")

figure(6)
title("T(T_{in})")

%% Zmiana 2. zakłócenia:
labels = [];

figure(7)
figure(8)
for i=1:7
	
	c0 = [C_Ain; F_C];
	d0 = [T_in; T_Cin*coef_vals(i)];
	
	[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);
	
	figure(7)
	hold on
	plot(t, x(:, 1), 'Color', color_formats(i, :))
	figure(8)
	hold on
	plot(t, x(:, 2), 'Color', color_formats(i, :))
	
	labels = [labels, "Model nieliniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end

lin_x0 = x0 - [C_A; T];

for i=1:7
	c0 = [C_Ain; F_C];
	d0 = [T_in; T_Cin*coef_vals(i)];
	lin_c0 = c0 - [C_Ain; F_C];
	lin_d0 = d0 - [T_in; T_Cin];
	
	[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
	lin_x(:, 1) = lin_x(:, 1) + C_A;
	lin_x(:, 2) = lin_x(:, 2) + T;
	
	figure(7)
	hold on
	plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(8)
	hold on
	plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
	
	labels = [labels, "Model liniowy, F_C zmnienione " + coef_vals(i) + " razy"];
end
figure(7)
% legend(labels)
title("C_A(T_{Cin})")

figure(8)
title("T(T_{Cin})")

% figure(9)
% step(c2d(tf(ss(A,[B E], eye(2), zeros(2,4))),1))

function dx = reactor(t, x, input, disturbance)
	global Ro Ro_c c_p c_pc k_0 E_R h a b;
	global V F_in F;

	d_Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1));
	d_T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)));

	dx = [d_Ca; d_T];
end

function dx = lin_reactor(t, x, input, disturbance)
	global A B E;
% 	global Ro Ro_c c_p c_pc k_0 E_R h a b;
% 	global V F_in F C_Ain F_C C_A T T_in T_Cin;
% 	d_Ca = (F_in*C_Ain/V - F*C_A/V - k_0*exp(-E_R/T)*C_A)/60;
% 	d_T = (F_in*T_in/V - F*T/V + h*k_0*exp(-E_R/T)*C_A/(Ro*c_p) - a*(F_C)^(b+1)/((F_C+ a*(F_C)^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(T-T_Cin))/60;

% 	dx = [d_Ca;d_T] + A*x + B*input;
	dx = A*x + B*input + E*disturbance;
end