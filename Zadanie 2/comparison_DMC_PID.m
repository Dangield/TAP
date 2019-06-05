clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

obj = NonlinearReactor();

load("data/setPoints.mat")
workpoint = struct('u', [2; 15], 'y', [0.2646; 393.9521], 'd', [323; 365]);
obj.resetToWorkPoint(workpoint);

umin = [0.1; 0.1];
umax = [4; 30];
dumax = [1; 2];

D = 800;
N = 100;
Nu = 80;
lambda = [2 0.1];
psii = 1;
sim_length = 6000;


load('./data/s.mat', 's');
dmc = DMC_Regulator(obj, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

params = {struct('Kp', 2, 'Ti', 0.03, 'Td', 2), struct('Kp', 0.1, 'Ti', 0.5, 'Td', 0.1)};
pid = PID_Regulator(obj, params, []);

u_dmc = workpoint.u.*ones(obj.nu, sim_length);
y_dmc = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y_dmc(:, k) = output;
    control = dmc.calculate(output, setPoints(:, k));
    u_dmc(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

obj.resetToWorkPoint(workpoint);

u_pid = workpoint.u.*ones(obj.nu, sim_length);
y_pid = workpoint.y.*ones(obj.ny, sim_length);

for k = 1:sim_length
    output = obj.getOutput();
    y_pid(:, k) = output;
    control = flipud(pid.calculate(output, setPoints(:, k)));
    u_pid(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

titles = ["C_A", "T"];
figure;
for i = 1:2
    subplot(2, 1, i);
        hold on;
        stairs(setPoints(i, :), 'b');
		stairs(y_dmc(i, :), 'r');
		stairs(y_pid(i, :), 'm');
		title(titles(i));
		legend("Wartosci zadane", "DMC Analityczny", "PID")
end

titles = ["C_Ain", "F_c"];
figure;
for i = 1:2
    subplot(2, 1, i);
		hold on
		stairs(u_dmc(i, :), 'r');
		stairs(u_pid(i, :), 'm');
		legend("DMC Analityczny", "PID")
end

