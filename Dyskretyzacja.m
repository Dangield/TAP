

c0 = [C_Ain*1.1; F_C];
d0 = [T_in; T_Cin];
lin_c0 = c0 - [C_Ain; F_C];
lin_d0 = d0 - [T_in; T_Cin];

%%% Algorytm Rungego-Kutty
for Tp=[0.1, 0.5, 1, 2, 3]
    d_T = [0];
    d_C_A = [0];
    d_Y = [d_T; d_C_A];
    T_s = [];
    d_Y_l = [];
    
    for i=0:Tp:15
        k1 = Tp*lin_reactor(0, d_Y, c0, d0);
        k2 = Tp*lin_reactor(0, d_Y + k1/2, c0 + [Tp; Tp]/2, d0 + [Tp; Tp]/2);
        k3 = Tp*lin_reactor(0, d_Y + k2/2, c0 + [Tp; Tp]/2, d0 + [Tp; Tp]/2);
        k4 = Tp*lin_reactor(0, d_Y + k3, c0 + [Tp; Tp], d0 + [Tp; Tp]);
        
        d_Y = (k1 + 2*k2 + 2*k3 + k4)/6;
        d_Y_l = [d_Y_l, d_Y];
        T_s = [T_s; i];
    end
    
    disp("=================================")
    Tp
    size(d_Y_l)
    size(T_s)
    
    figure(1)
    hold on
    plot(T_s, d_Y_l(1, :)')
    
    figure(2)
    hold on
    plot(T_s, d_Y_l(2, :)')
end
[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [0 15], [0, 0]);

figure(1)
hold on
plot(lin_t, lin_x(:, 1))
figure(2)
hold on
plot(lin_t, lin_x(:, 2))

function dx = lin_reactor(t, x, input, disturbance)
	global A B E;
% 	global Ro Ro_c c_p c_pc k_0 E_R h a b;
% 	global V F_in F C_Ain F_C C_A T T_in T_Cin;
% 	d_Ca = (F_in*C_Ain/V - F*C_A/V - k_0*exp(-E_R/T)*C_A)/60;
% 	d_T = (F_in*T_in/V - F*T/V + h*k_0*exp(-E_R/T)*C_A/(Ro*c_p) - a*(F_C)^(b+1)/((F_C+ a*(F_C)^b/(2*Ro_c*c_pc))*(V*Ro*c_p))*(T-T_Cin))/60;

% 	dx = [d_Ca;d_T] + A*x + B*input;
	dx = A*x + B*input + E*disturbance;
end