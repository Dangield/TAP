classdef Reactor < AbstractObject
	properties
		u, y, d
		uk, yk, dk
	end
	
	methods
		function self = Reactor()
			ny = 2; nu = 2; nd = 2; Ts = 1/60;
			self@AbstractObject(ny, nu, nd, Ts);
		end
		
		function output = getOutput(self)
			output = self.yk;
		end
		
		function setControl(self, control)
			self.uk = control;
		end
			
		function nextIteration(self)
			self.shiftArrays();
			self.simulate();
		end
		
		function shiftArrays(self)
			self.u = circshift(self.u, [0 1]);
			self.y = circshift(self.y, [0 1]);
            self.d = circshift(self.y, [0 1]);
			
			self.u(:, 1) = self.uk;
			self.y(:, 1) = self.yk;
            self.d(:, 1) = self.dk;
		end
		
		function simulate(self)
			global Ro Ro_c c_p c_pc k_0 E_R h a b;
			global V F_in F;
			input = self.uk;
			x = self.y(:, 1);
			disturbance = self.d(:, 1);

			Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))*self.Ts + x(1);
			T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))*self.Ts + x(2);
			self.yk = [Ca; T];
		end
		
		function resetToWorkPoint(self, workPoint)
			self.uk = workPoint.u;
			self.yk = workPoint.y; 
			self.dk = workPoint.d; 

			self.u = self.uk*ones(1, 5);
			self.y = self.yk*ones(1, 5);
			self.d = self.dk*ones(1, 5);
		end
	end
end

