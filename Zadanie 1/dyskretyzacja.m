close all
clear
addpath("../")

consts
linearyzacja;
clc

t0 = 0;
tfinal = 10;
x0 = [C_A; T];
d0 = [T_in; T_Cin];
c0 = [C_Ain * 1.2 C_Ain; F_C, F_C * 1.4];

outputs = ["C_{A}", "T"];
inputs = ["C_{Ain}", "F_C"];

linearReactor = LinearReactor();

for Ts  = [0.1 0.05 0.01]
	figure
	for inputJump = 1:2
		[t, x] = linearReactor.simulate(x0, c0(:, inputJump), d0, t0, tfinal);

		discreteTransferFunction = c2d(transferFunction, Ts);
		discreteSpaceState       = ss(discreteTransferFunction);
		disctreteReactor         = DiscreteReactor(discreteSpaceState);
		[t_d, x_d] = disctreteReactor.simulate(x0, c0(:, inputJump), d0, t0, tfinal);

		for output = 1:2
			subplot(2, 2, (output-1)*2+inputJump)
				hold on
				plot(t, x(:,output))
				stairs(t_d, x_d(:, output))
				title(outputs(output) + ", skok na " + inputs(inputJump) + " = " + c0(inputJump, inputJump));
				xlabel('t[min]')
				ylabel(outputs(output))
				legend("Model ci¹g³y liniowy", "Model dyskretny, T_s = " + Ts, 'Location', 'best');
		end
	end
end