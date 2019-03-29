close all
clear
clc
addpath("../")

consts

global A B E;
A = [
	-F/V - k_0*exp(-E_R/T), -k_0*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A;
	h*k_0*exp(-E_R/T)/Ro/c_p, -F/V+(k_0*h*V*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A -a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc)))/(V*Ro*c_p)
	];

B = [
	F_in/V, 0;
	0, -((F_C^b*a*(T - T_Cin)*(b + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro))...
	- (F_C^(b + 1)*a*(T - T_Cin)*((F_C^(b - 1)*a*b)/(2*c_pc*Ro_c) + 1))/(F_C + (F_C^b*a)/(2*c_pc*Ro_c))^2)/(V*c_p*Ro)
	];

E = [
	0, 0;
	F_in/V, a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
	];

stateSpaceModel  = ss(A,[B E], eye(2), zeros(2,4))
transferFunction = tf(stateSpaceModel)