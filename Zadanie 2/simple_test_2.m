clear;
% close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

obj = SimpleObject2();

workpoint = struct('u', [0; 0], 'y', [0; 0]);
obj.resetToWorkPoint(workpoint);

umin = [-100; -100];
umax = [100; 100];
dumax = [100; 100];

D = 400;
N = 100;
Nu = 100;
lambda = [0 0];
psii = 1;
sim_length = 100;

load('./data/s2.mat', 's');
reg = Numeric_DMC_Regulator(obj, workpoint, s, D, N, Nu, lambda, psii, umin, umax, dumax);

u = workpoint.u.*ones(obj.nu, sim_length);
y = workpoint.y.*ones(obj.ny, sim_length);

setPoints = y;
setPoints(:, 50:end) = 1;

for k = 1:sim_length
    output = obj.getOutput();
    y(:, k) = output;
    control = reg.calculate(output, setPoints(:, k));
    u(:, k) = control';
    obj.setControl(control);
    obj.nextIteration();
end

% figure;
%     hold on
% 	stairs(u', 'r');
% 	title("u")
% 	
% figure;
% 	hold on;
% 	stairs(setPoints', 'b');
% 	stairs(y', 'r');
% 	title("y")
% 	legend("yzad", "y")
	
figure;
for i = 1:2
    subplot(2, 1, i);
		stairs(u(i, :), 'r');
end

figure;
for i = 1:2
    subplot(2, 1, i);
        hold on;
        stairs(setPoints(i, :), 'b');
		stairs(y(i, :), 'r');
end
