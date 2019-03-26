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
C_Ain = 2;
F_C = 15;
T_in = 323;
T_Cin = 365;
C_A = 0.26;
T = 393.9;

t0 = 0
tfinal = 100
% x0 = [400; 450]
% x0 = [260; 393.9]
x0 = [0.26; 393.9]

[t, x] = ode23(@reactor, [t0 tfinal], x0);

plot(t, x)
title("Values")
legend("C_A", "T")

R_U = (F_in*C_Ain - F*C_A - V*k_0*exp(-E_R/T)*C_A)/V

R_D = (F_in*Ro*c_p*T_in - F*Ro*c_p*T + V*h*k_0*exp(-E_R/T)*C_A - a*(F_C)^(b+1)*(T-T_Cin)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc)))/(V*Ro*c_p)

% [F_in*Ro*c_p*T_in; - F*Ro*c_p*T; V*h*k_0*exp(-E_R/T)*C_A; - a*(F_C)^(b+1)*(T-T_Cin)/(F_C+a*(F_C)^b/(2*Ro_c*c_pc))]


function dx = reactor(t, x)

    global Ro Ro_c c_p c_pc k_0 E_R h a b;
    global V F_in F C_Ain F_C T_in T_Cin;

    
    d_Ca = F_in*C_Ain/V - F*x(1)/V - k_0*exp(-E_R/x(2))*x(1);
    d_T = F_in*T_in/V - F*x(2)/V + h*k_0*exp(-E_R/x(2))*x(1)/(Ro*c_p) - a*(F_C)^(b+1)/((F_C + a*(F_C)^b/(2*Ro_c+c_pc))*(V*Ro*c_p))*(x(2)-T_Cin);
    
%     if t < 10000
%         t
%         C_A_t = x(1)
%         T_t = x(2)
%         disp('============================================================================')
%     end
    
    dx = [d_Ca; d_T];
end