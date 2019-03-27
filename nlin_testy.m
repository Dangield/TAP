

dC_Anlin = [];
dC_Alin = [];

dT_nlin = [];
dT_lin = [];

C_A_vals = [];
for C_A_test=0.1:0.01:0.4
    val = reactor([C_A_test; T], [C_Ain; F_C], [T_in; T_Cin]);
    lin_val = lin_reactor([C_A_test-C_A; 0], [0; 0], [0; 0]);
	dC_Anlin = [dC_Anlin; val(1)];
	dC_Alin = [dC_Alin; (F_in*C_Ain/V - F*C_A/V - k_0*exp(-E_R/T)*C_A/V)/60 + lin_val(1)];
    
	dT_nlin = [dT_nlin; val(2)];
	dT_lin = [dT_lin; (F_in*C_Ain/V - F*C_A/V - k_0*exp(-E_R/T)*C_A/V) + lin_val(2)];
    
    C_A_vals = [C_A_vals; C_A_test];
end

subplot(2, 2, 1)
plot(C_A_vals, dC_Anlin)
hold on
plot(C_A_vals, dC_Alin)
title("d(dC_A/dt)/dC_A")

subplot(2, 2, 2)
plot(C_A_vals, dT_nlin)
hold on
plot(C_A_vals, dT_lin)
title("d(dT/dt)/dC_A")

dC_Anlin = [];
dC_Alin = [];

dT_nlin = [];
dT_lin = [];

dT_pp = reactor([C_A; T], [C_Ain; F_C], [T_in; T_Cin])
dT_pp = dT_pp(2)

T_vals = [];
for T_test=340:1:440
    val = reactor([C_A; T_test], [C_Ain; F_C], [T_in; T_Cin]);
    lin_val = lin_reactor([0; T_test-T], [0; 0], [0; 0]);
	dC_Anlin = [dC_Anlin; val(1)];
	dC_Alin = [dC_Alin; dT_pp + lin_val(1)];
    
	dT_nlin = [dT_nlin; val(2)];
	dT_lin = [dT_lin; dT_pp + lin_val(2)];
    
    T_vals = [T_vals; T_test];
end

subplot(2, 2, 3)
plot(T_vals, dC_Anlin)
hold on
plot(T_vals, dC_Alin)
title("d(dC_A/dt)/dT")

subplot(2, 2, 4)
plot(T_vals, dT_nlin)
hold on
plot(T_vals, dT_lin)
title("d(dT/dt)/dT")

function dx = reactor(x, input, disturbance)
	global Ro Ro_c c_p c_pc k_0 E_R h a b;
	global V F_in F;

	d_Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))/60;
	d_T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))/60;

	dx = [d_Ca; d_T];
end

function dx = lin_reactor(x, input, disturbance)
	global A B E;
% 	global Ro Ro_c c_p c_pc k_0 E_R h a b;
% 	global V F_in F C_Ain F_C C_A T T_in T_Cin;
% 	d_Ca = (F_in*C_Ain/V - F*C_A/V - k_0*exp(-E_R/T)*C_A)/60;
% 	d_T = (F_in*T_in/V - F*T/V + h*k_0*exp(-E_R/T)*C_A/(Ro*c_p) - a*(F_C)^(b+1)/((F_C+ a*(F_C)^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(T-T_Cin))/60;

% 	dx = [d_Ca;d_T] + A*x + B*input;
	dx = A*x + B*input + E*disturbance;
end