clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

obj = NonlinearReactor();

workpoint = struct('u', [2; 15], 'y', [0.2646; 393.9521], 'd', [323; 365]);
obj.resetToWorkPoint(workpoint);

sim_length = 1000;
u = (workpoint.u.*ones(obj.nu, sim_length));
y = (workpoint.y.*ones(obj.ny, sim_length));

u(1, 500:1000) = 2.4;

for k = 1:sim_length
    y(:, k) = obj.getOutput();
    obj.setControl(u(:, k));
    obj.nextIteration();
end

figure;
subplot(2, 2, 1);
    stairs(u(1, :));
    title('U1');
subplot(2, 2, 2);
    stairs(u(2, :));
    title('U2');
subplot(2, 2, 3);
    stairs(y(1, :));
    legend('Y1');
subplot(2, 2, 4);
    stairs(y(2, :));
    legend('Y2');