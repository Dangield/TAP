global Ro Ro_c c_p c_pc k_0 E_R h a b;
global V F_in F C_Ain F_C T_in T_Cin;

% % Stałe
% Ro = 1000;
% Ro_c = 1000;
% c_p = 1000;
% c_pc = 1000;
% k_0 = 10^10/60;
% E_R = 8330.1;
% h = 130*10^3;
% a = 1.678*10^6;
% b = 0.5;
% 
% % Wartości sterujące
% V = 1;
% F_in = 1/60;
% F = 1/60;
% C_Ain = 2000;
% F_C = 15/60;
% T_in = 323;
% T_Cin = 365;
% C_A = 260;
% T = 393.9;


%Stałe
Ro = 1e6;
Ro_c = 1e6;
c_p = 1;
c_pc = 1;
k_0 = 1e10;
E_R = 8330.1;
h = 130e6;
a = 1.678e6;
b = 0.5;

%Wartości sterujące
V = 1;
F_in = 1;
F = 1;
% Sterowanie
C_Ain = 2;
F_C = 15;
% Zakłócenie
T_in = 323;
T_Cin = 365;
C_A = 0.26;
T = 393.9;

% linearyzacja:
global A B E
A = [
    -F/V - k_0*exp(-E_R/T), -k_0*(-E_R)*(-1/T^2)*exp(-E_R/T)*C_A;
    h*k_0*exp(-E_R/T)/Ro/c_p, -F/V-a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
    ]/60;
eig(A)
B = [-F/V, 0;
    0, -a*(b+1)*(F_C)^b*(1+a*(F_C)^(b-1)/(2*Ro_c*c_pc))/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))^2/(V*Ro*c_p)
    ]/60;
E = [0, 0;
    F_in/V, -a*(F_C)^(b+1)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))/(V*Ro*c_p)
    ]/60;



t0 = 0;
tfinal = 500;
x0 = [0.23; 393.9];
c0 = [2; 15];
d0 = [323; 365];

% x0 = [400; 450]
% x0 = [260; 393.9]
% x0 = [0.26; 393.9]

[t, x] = ode45(@(t, x) reactor(t, x, c0, d0), [t0 tfinal], x0);

lin_x0 = x0 - [C_A; T];
lin_c0 = c0 - [C_Ain; F_C];
lin_d0 = d0 - [T_in; T_Cin];
[lin_t, lin_x] = ode45(@(t, x) lin_reactor(t, x, lin_c0, lin_d0), [t0 tfinal], lin_x0);
lin_x(:, 1) = lin_x(:, 1) + C_A;
lin_x(:, 2) = lin_x(:, 2) + T;

plot(t, x(:, 1))
hold on
plot(lin_t, lin_x(:, 1))
hold off
title("C_A")
figure()
plot(t, x(:, 2))
hold on
plot(lin_t, lin_x(:, 2))
hold off
title("T")

function dx = reactor(t, x, input, disturbance)

    global Ro Ro_c c_p c_pc k_0 E_R h a b;
    global V F_in F;

    
%     d_Ca = (F_in*C_Ain/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))/60;
%     d_T = (F_in*T_in/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(F_C)^(b+1)/((F_C + a*(F_C)^b/(2*Ro_c+c_pc))*(V*Ro*c_p))*(x(2)-T_Cin))/60;
    
    d_Ca = (F_in*input(1)/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1))/60;
    d_T = (F_in*disturbance(1)/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(input(2))^(b+1)/((input(2)+ a*(input(2))^b/(2*Ro_c+c_pc))*(V*Ro*c_p))*(x(2)-disturbance(2)))/60;
    
%     if t < 10000
%         t
%         C_A_t = x(1)
%         T_t = x(2)
%         disp('============================================================================')
%     end
    
    dx = [d_Ca; d_T];
end

function dx = lin_reactor(t, x, input, disturbance)

    global A B E;
    
    dx = A*x + B*input + E*disturbance;
end