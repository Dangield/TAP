clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

obj = SimpleObject();

workpoint = struct('u', 0, 'y', 0);
obj.resetToWorkPoint(workpoint);

umin = -100;
umax = 100;
dumax = 100;

D = 6;
N = 100;
Nu = 10;
lambda = 0;
psii = 1;
sim_length = 100;

load('./data/simpleS.mat', 's');
reg = DMC_Regulator(obj, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setPoints = y;
setPoints(50:end) = 1;

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = reg.calculate(output, setPoints(:, k));
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

figure;
	stairs(u(1, :), 'r');
	title("u")
	
figure;
	hold on;
	stairs(setPoints(1, :), 'b');
	stairs(y(1, :), 'r');
	title("y")
	legend("yzad", "y")
	
