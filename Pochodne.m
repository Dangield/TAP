syms F_in F_c F V k_0 E R 
syms ro ro_c c_p c_pc T_in T_Cin T h a b E_R
syms C_A T C_Ain F_c Ca0 T0

dCadt=(F_in*C_Ain -F*C_A - V*k_0*exp(-E_R/T)*C_A)/V %(V przerzucamy na drug¹ stronê
dTdt= (F_in*ro*c_p*T_in - F*ro*c_p*T +...
    V*h*k_0*exp(-E_R/T)*C_A -((F_c^(b + 1)*a*(T-T_Cin))/(F_c + ((F_c^b*a)/(2*c_pc*ro)))))/(V*c_p*ro) %%V*c_p*ro przerzucamy na drug¹ stronê


%Wyliczamy pochodne cz¹stkowe
disp('dCadt po C_A');
dCadt_Ca=diff(dCadt,C_A) 
disp('dCadt po T');
dCadt_T=diff(dCadt,T) 
disp('dCadt po C_Ain');
dCadt_C_Ain=diff(dCadt,C_Ain)
disp('dCadt po Fc');
dCadt_Fc=diff(dCadt,F_c) 

%---------------------------------
disp('dTdt po C_A');
dTdt_Ca=diff(dTdt,C_A)
disp('dTdt po T');
dTdt_T=diff(dTdt,T)
disp('dTdt po C_Ain');
dTdt_C_Ain=diff(dTdt,C_Ain)
disp('dTdt po Fc');
dTdt_Fc=diff(dTdt,F_c)


