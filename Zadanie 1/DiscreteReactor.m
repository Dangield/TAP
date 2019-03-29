classdef DiscreteReactor

	properties
		A;
		B;
		C;
		D;
		Ts;
	end
	
	methods
		function obj = DiscreteReactor(discreteSpaceState)
			obj.A  = discreteSpaceState.A;
			obj.B  = discreteSpaceState.B;
			obj.C  = discreteSpaceState.C;
			obj.D  = discreteSpaceState.D;
			obj.Ts = discreteSpaceState.Ts;
		end
		
		function [t, y] = simulate(obj, x0, c0, d0, t0, tfinal)
			global C_Ain F_C T_in T_Cin C_A T;
			lin_x0 = x0 - [C_A; T];
			lin_c0 = c0 - [C_Ain; F_C];
			lin_d0 = d0 - [T_in; T_Cin];
			
			x = zeros(size(obj.A, 1), 1);
			y = lin_x0;
			t = t0:obj.Ts:tfinal;
			for i = 2:size(t, 2)
				x(:, i) = obj.A*x(:, i-1) + obj.B*[lin_c0; lin_d0];
				y(:, i) = obj.C*x(:, i);
			end
			y = y';
			y(:, 1) = y(:, 1) + C_A;
			y(:, 2) = y(:, 2) + T;
		end
	end
end

