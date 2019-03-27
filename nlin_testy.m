


C_Anlin = [];
C_Alin = [];
for C_Ain_=0.1:0.1:5
	C_Anlin = [C_Anlin;F_in*C_Ain_/(F+ V*k_0*exp(-E_R/T))];
	C_Alin = [C_Alin; C_A+60*B(1, 1)*(C_Ain_-C_Ain)];
end

plot(C_Anlin)
hold on
plot(C_Alin)

	
	