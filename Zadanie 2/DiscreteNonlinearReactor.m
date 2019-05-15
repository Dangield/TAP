classdef DiscreteNonlinearReactor < handle

	properties
		Ts;
		u;
		y;
		d;
		t;
		Tp = 1/60;
		currentU;
	end
	
	methods
		function obj = DiscreteNonlinearReactor()
			global C_Ain F_C  C_A T  T_in T_Cin;
			obj.u = [C_Ain, F_C];
			obj.currentU = [C_Ain, F_C];
			obj.y = [C_A, T];
			obj.d = [T_in, T_Cin];
			obj.t = 0;
		end
		
		
		function simulate(obj)
			global Ro Ro_c c_p c_pc k_0 E_R h a b;
			global V F_in F;
			input = obj.currentU;
			x = obj.y(end, :);
			disturbance = obj.d(end, :);

			Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))*obj.Tp + x(1);
			T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))*obj.Tp + x(2);
			obj.u = [obj.u; obj.currentU];
			obj.y = [obj.y; Ca T];
			obj.t = [obj.t; obj.t(end) + obj.Tp];
		end
	end
end
