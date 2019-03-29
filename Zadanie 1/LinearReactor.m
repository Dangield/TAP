classdef LinearReactor
	properties
		
	end
	
	methods
		function obj = LinearReactor()
		end
		
		function [t, lin_x] = simulate(obj, x0, c0, d0, t0, tfinal)
			global C_Ain F_C T_in T_Cin C_A T;
			lin_x0 = x0 - [C_A; T];
			lin_c0 = c0 - [C_Ain; F_C];
			lin_d0 = d0 - [T_in; T_Cin];
			[t, x] = ode45(@(t, x) obj.differential(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
			lin_x(:, 1) = x(:, 1) + C_A;
			lin_x(:, 2) = x(:, 2) + T;
		end
		
		function dx = differential(obj, t, x, input, disturbance)
			global A B E;
			dx = A*x + B*input + E*disturbance;
		end
	end
end

