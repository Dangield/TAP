clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

% TODO - use LINEAR REACTOR !!!
obj = SimpleObject();

workpoint = struct('u', 0, 'y', 0);
obj.resetToWorkPoint(workpoint);

sim_length = 100;

s = cell(obj.ny, obj.nu);

start = 10;

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

save('./data/simpleS.mat', 's');
figure;
    stairs(s{1, 1});