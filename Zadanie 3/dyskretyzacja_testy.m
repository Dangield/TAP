close all
clear all
clc
addpath("../")

consts

C_A_w = 0.26;
T_w = 393.9;

C_Ain_w = 2;
F_C_w= 15;

T_in_w = 323;
T_Cin_w = 365;

Tp = 1/60;

last_x = [C_A_w; T_w];
total_x = last_x';
t = [0];
for i=0:Tp:50
    last_x = reactor(last_x, [C_Ain_w, F_C_w], [T_in_w, T_Cin_w], Tp);
    total_x = [total_x; last_x'];
    t = [t; t(end) + Tp];
end

figure(1)
stairs(t, total_x(:, 1))
title('C_A')
figure(2)
stairs(t, total_x(:, 2))
title('T')
total_x(end, :)

function xn = reactor(x, input, disturbance, Tp)
	global Ro Ro_c c_p c_pc k_0 E_R h a b;
	global V F_in F;

	Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))*Tp + x(1);
	T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))*Tp + x(2);

	xn = [Ca; T];
end