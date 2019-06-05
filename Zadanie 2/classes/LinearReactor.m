classdef LinearReactor < AbstractObject
	properties
		u, y, d
		uk, yk, dk
		
		A;
		B;
		C;
		D;
	end
	
	methods
		function self = LinearReactor()
			ny = 2; nu = 2; nd = 2; Ts = 1/60;
			self@AbstractObject(ny, nu, nd, Ts);
			
			load("data/discreteLinearParameters.mat");
			self.A  = A;
			self.B  = B;
			self.C  = C;
			self.D  = D;
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
            self.d = circshift(self.d, [0 1]);
			
			self.u(:, 1) = self.uk;
			self.y(:, 1) = self.yk;
            self.d(:, 1) = self.dk;
		end
		
		function simulate(self)
			global C_Ain F_C T_in T_Cin C_A T;
			lin_x0 = self.yk - [C_A; T];
			lin_c0 = self.uk' - [C_Ain; F_C];
			lin_d0 = self.dk - [T_in; T_Cin];
			
			x = self.A * lin_x0 + self.B*[lin_c0; lin_d0];
			y = self.C * x;
            
			y(1) = y(1) + C_A;
			y(2) = y(2) + T;
			
			self.yk = y;
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

