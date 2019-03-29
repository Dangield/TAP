classdef Reactor
	properties
		
	end
	
	methods
		function obj = Reactor()
		end
		
		function [t, x] = simulate(obj, x0, c0, d0, t0, tfinal)
			[t, x] = ode45(@(t, x) obj.differential(t, x, c0, d0), [t0 tfinal], x0);
		end
		
		function dx = differential(obj, t, x, input, disturbance)
			global Ro Ro_c c_p c_pc k_0 E_R h a b;
			global V F_in F;

			d_Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1));
			d_T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)));

			dx = [d_Ca; d_T];
		end
	end
end

