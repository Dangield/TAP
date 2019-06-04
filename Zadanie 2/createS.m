clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

% TODO - use LINEAR REACTOR !!!
obj = LinearReactor();

workpoint = struct('u', [2; 15], 'y', [0.2646; 393.9521], 'd', [323; 365]);
obj.resetToWorkPoint(workpoint);

sim_length = 1000;

s = cell(obj.ny, obj.nu);

start = 100;

for n = 1:obj.nu
    u = workpoint.u.*ones(obj.nu, sim_length);
    y = workpoint.y.*ones(obj.ny, sim_length);
    u(n, start:end) = workpoint.u(n) + 1;
    for k = 1:sim_length
        y(:, k) = obj.getOutput();
        obj.setControl(u(:, k));
        obj.nextIteration();
    end
    obj.resetToWorkPoint(workpoint);
    for m = 1:obj.ny
        s{m, n} = y(m, start+1:end) - y(m, start);
    end
end

save('./data/s.mat', 's');
figure;
subplot(2, 2, 1);
    hold on;
    stairs(s{1, 1});
    title('U1');
subplot(2, 2, 2);
    stairs(s{2, 1});
    title('U1');
	
subplot(2, 2, 3);
    stairs(s{1, 2});
    title('U2');
subplot(2, 2, 4);
    stairs(s{2, 2});
    title('U2');