close all
clear
addpath("../")

consts
linearyzacja
clc

% Symulacja
t0 = 0;
tfinal = 20;
x0 = [C_A; T];
u0 = [C_Ain; F_C; T_in; T_Cin];
coef_vals  = [0.4, 0.6, 0.8, 1, 1.2, 1.4, 1.6];

reactor       = Reactor();
linearReactor = LinearReactor();

simResults = cell(4, 2);

for inputIterator = 1:4
	u = u0;
	for i=1:size(coef_vals, 2)
		u(inputIterator) = u0(inputIterator) * coef_vals(i);
		c0 = u(1:2);
		d0 = u(3:4);
		
		[t, x] = reactor.simulate(x0, c0, d0, t0, tfinal);
		[tLin, xLin] = linearReactor.simulate(x0, c0, d0, t0, tfinal);
		
		simResults{inputIterator, 1}{i} = [t, x];
		simResults{inputIterator, 2}{i} = [tLin, xLin];
	end
end

% Rysowanie
color_formats = [[1 0 0]; [0 1 0]; [0 0 1]; [0 0.7 0.7]; [0.7 0 0.7]; [0.7 0.7 0]; [0 0 0]];
inputLabels  = ["C_{Ain}", "F_C", "T_{in}", "T_{Cin}"];
outputLabels = ["C_A", "T"];

for inputIterator = 1:4
	inputLabel = inputLabels(inputIterator);
	for outputIterator = 1:2
		legendLabels = [];
		figure
			hold on
			title(outputLabels(outputIterator) + "(" + inputLabel + ")")
			for jumpIterator=1:size(coef_vals, 2)
				t = simResults{inputIterator, 1}{jumpIterator}(:, 1);
				x = simResults{inputIterator, 1}{jumpIterator}(:, 1 + outputIterator);
				t_lin = simResults{inputIterator, 2}{jumpIterator}(:, 1);
				x_lin = simResults{inputIterator, 2}{jumpIterator}(:, 1 + outputIterator);
				
				plot(t, x, 'Color', color_formats(jumpIterator, :))
				plot(t_lin, x_lin, 'Color',  color_formats(jumpIterator, :), 'LineStyle', '--')

				legendLabels = [legendLabels, "Model nieliniowy, " + inputLabel + " = " + u0(inputIterator) * coef_vals(jumpIterator), ...
					"Model liniowy, " + inputLabel + " = " + u0(inputIterator) * coef_vals(jumpIterator)];
			end
			legend(legendLabels,'Location','eastoutside');
	end
end