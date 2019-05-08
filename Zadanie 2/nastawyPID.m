% close all
clear all
clc
addpath("../")

consts

Tp = 1/60;

K = 0;          % Gain
Ti = 0;         % Integration time
Kd = 0;         % Derivation gain
Td = 0;         % Derivation time
Tp = 1/60;      % Sample time
Dir = 0;        % 1 - direct (SP-PV), 0 - indirect (PV-SP)
Hlim = 100;     % high limit
Llim = -100;    % low limit
AutoMan = 1;    % MAN when AutoMan equals to 0
ManVal = 1;     % output in MAN mode;

r = DiscreteNonlinearReactor(Tp);
p1 = classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, AutoMan, ManVal);
p2 = classPID(K, Ti, Kd, Td, Tp, Hlim, Llim, AutoMan, ManVal);

%SetAutoMan(obj, AutoMan, ManVal)
p1.SetAutoMan(1,1)
p2.SetAutoMan(1,1)

%PID
%reTune(K, Ti, Kd, Td)
p1.reTune(15, 0.01, 40, 0.1)
p2.reTune(0.01, 3, 0.2, 0.01)

global C_A T C_Ain F_C;

yzad = [C_A, T];
u = [C_Ain F_C];

yzad_log = yzad;


for i=0:Tp:250
	i
	u2 = p1.calc(r.y(end, 1), yzad(1)) + F_C;
	u1 = p2.calc(r.y(end, 2), yzad(2)) + C_Ain;
	if (i >= 50)
		yzad(2) = 400;
	end
	
	if (i >= 150)
		yzad(1) = 0.2;
	end
	r.currentU = [u1 u2];
    r.simulate()
	yzad_log = [yzad_log; yzad];
end

figure()
	subplot(2, 1, 1)
		hold on
		stairs(r.t, r.y(:, 1))
		stairs(r.t, yzad_log(:, 1))
		title('C_A')
	subplot(2, 1, 2)
		hold on
		stairs(r.t, r.y(:, 2))
		stairs(r.t, yzad_log(:, 2))
		title('T')