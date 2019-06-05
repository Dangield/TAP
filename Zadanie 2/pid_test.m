clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

sim_length = 1000;

consts


% stateSpaceModel  = ss(A,[B E], eye(2), zeros(2,4))
% noPerturbsStateSpaceModel = ss(A, B, eye(2), zeros(2,2));
% transferFunction = tf(noPerturbsStateSpaceModel)


params = {struct('Kp', 2, 'Ti', 0.03, 'Td', 2), struct('Kp', 0.1, 'Ti', 0.5, 'Td', 0.1)};

%% Ma≈Ça zmiana w C_A
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [0.29; 393.9521].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zadane C_A=0.29, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zadane C_A=0.29, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zadane C_A=0.29, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zadane C_A=0.29, model liniowy')

%% Du≈ºa zmiana w C_A
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [0.37; 393.9521].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zadane C_A=0.37, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zadane C_A=0.37, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zadane C_A=0.37, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zadane C_A=0.37, model liniowy')

%% Ma≈Ça zmiana w T
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; 395].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zadane T=395, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zadane T=395, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zadane T=395, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zadane T=395, model liniowy')

%% Du≈ºa zmiana w T
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; 425].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zadane T=425, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zadane T=425, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zadane T=425, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zadane T=425, model liniowy')

%% Ma≈Ça zmiana w T_{in}
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in + 2; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; T].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zak≥Ûcenie T_{in}=325, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zak≥Ûcenie T_{in}=325, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zak≥Ûcenie T_{in}=325, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zak≥Ûcenie T_{in}=325, model liniowy')

%% Du≈ºa zmiana w T_{in}
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in + 20; T_Cin]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; T].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zak≥Ûcenie T_{in}=343, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zak≥Ûcenie T_{in}=343, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zak≥Ûcenie T_{in}=343, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zak≥Ûcenie T_{in}=343, model liniowy')

%% Ma≈Ça zmiana w T_Cin
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin + 2]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; T].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zak≥Ûcenie T_Cin=367, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zak≥Ûcenie T_Cin=367, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zak≥Ûcenie T_Cin=367, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zak≥Ûcenie T_Cin=367, model liniowy')

%% Du≈ºa zmiana w T_Cin
% obj = NonlinearReactor();
obj = LinearReactor();

workpoint = struct('u', [C_Ain; F_C], 'y', [C_A; T], 'd', [T_in; T_Cin + 20]);
obj.resetToWorkPoint(workpoint);

reg = PID_Regulator(obj, params, []);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setY = [C_A; T].*ones(obj.ny, sim_length);
% setY = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = flipud(reg.calculate(output, setY(:, k), []) + [F_C; C_Ain]);
    if k < 100
    end
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure
subplot(2, 1, 1)
stairs(u(1, :), 'r');
title('C_{Ain}, Zak≥Ûcenie T_Cin=385, model liniowy')

subplot(2, 1, 2)
stairs(u(2, :), 'r');
title('F_{C}, Zak≥Ûcenie T_Cin=385, model liniowy')

figure;
subplot(2, 1, 1);
hold on;
stairs(setY(1, :), 'b');
stairs(y(1, :), 'r');
title('C_A, Zak≥Ûcenie T_Cin=385, model liniowy')

subplot(2, 1, 2);
hold on;
stairs(setY(2, :), 'b');
stairs(y(2, :), 'r');
title('T, Zak≥Ûcenie T_Cin=385, model liniowy')