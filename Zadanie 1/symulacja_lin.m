close all
clear
clc
addpath("../")

consts
linearyzacja

% Zmienne do plotowania
labels       = [];
linearLabels = [];
color_formats = [[1 0 0]; [0 1 0]; [0 0 1]; [0 0.7 0.7]; [0.7 0 0.7]; [0.7 0.7 0]; [0 0 0]];

% Zmienne do symulacji
t0 = 0;
tfinal = 20;
x0 = [C_A; T];
d0 = [T_in; T_Cin];
u0 = [C_Ain; F_C];
coef_vals  = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6];

reactor       = Reactor();
linearReactor = LinearReactor();

for i=1:size(coef_vals, 2)
	c0 = [C_Ain*coef_vals(i); F_C];
	
	[t, x]         = reactor.simulate(x0, c0, d0, t0, tfinal);
	[lin_t, lin_x] = linearReactor.simulate(x0, c0, d0, t0, tfinal);
	
	labels = [labels, "Model nieliniowy, C_{Ain} = " + C_Ain*coef_vals(i), ...
		"Model liniowy,	C_{Ain} = " + C_Ain*coef_vals(i)];
	
	figure(1)
		hold on
		plot(t, x(:, 1), 'Color', color_formats(i, :))
		plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(2)
		hold on
		plot(t, x(:, 2), 'Color', color_formats(i, :))
		plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
end

figure(1)
	title("C_A")
	legend(labels);

figure(2)
	title("T")
	legend(labels);

% zmiana drugiego sterowania
labels = [];
coef_vals = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6];

for i=1:size(coef_vals, 2)
	c0 = [C_Ain; F_C*coef_vals(i)];
	
	[t, x] = reactor.simulate(x0, c0, d0, t0, tfinal);
	[lin_t, lin_x] = linearReactor.simulate(x0, c0, d0, t0, tfinal);
	
	labels = [labels, "Model nieliniowy, F_C = " + F_C*coef_vals(i), ...
		"Model liniowy, F_C = " + F_C*coef_vals(i)];
	
	figure(3)
		hold on
		plot(t, x(:, 1), 'Color', color_formats(i, :))
		plot(lin_t, lin_x(:, 1), 'Color',  color_formats(i, :), 'LineStyle', '--')
	figure(4)
		hold on
		plot(t, x(:, 2), 'Color', color_formats(i, :))
		plot(lin_t, lin_x(:, 2), 'Color',  color_formats(i, :), 'LineStyle', '--')
end

figure(3)
	title("C_A")
	legend(labels)
figure(4)
	title("T")
	legend(labels)
