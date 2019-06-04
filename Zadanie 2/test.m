clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts

linear = LinearReactor();
nonLinear = NonlinearReactor();

workpoint = struct('u', [2; 15], 'y', [0.2646; 393.9521], 'd', [323; 365]);
linear.resetToWorkPoint(workpoint);
nonLinear.resetToWorkPoint(workpoint);

sim_length = 1000;
ul = (workpoint.u.*ones(linear.nu, sim_length));
yl = (workpoint.y.*ones(linear.ny, sim_length));

unl = (workpoint.u.*ones(linear.nu, sim_length));
ynl = (workpoint.y.*ones(linear.ny, sim_length));

ul(1, 500:1000) = 2.4;
unl(1, 500:1000) = 2.4;

for k = 1:sim_length
    yl(:, k) = linear.getOutput();
    linear.setControl(ul(:, k));
    linear.nextIteration();

	ynl(:, k) = nonLinear.getOutput();
    nonLinear.setControl(unl(:, k));
    nonLinear.nextIteration();
end

figure;
subplot(2, 2, 1);
	hold on
    stairs(ul(1, :));
	stairs(unl(1, :));
    title('U1');
subplot(2, 2, 2);
	hold on
    stairs(ul(2, :));
	stairs(unl(2, :));
    title('U2');
subplot(2, 2, 3);
	hold on
    stairs(yl(1, :));
	stairs(ynl(1, :));
    legend('Y1 - liniowe', 'Y1 - nieliniowe');
subplot(2, 2, 4);
	hold on
    stairs(yl(2, :));
	stairs(ynl(2, :));
    legend('Y2 - liniowe', 'Y2 - nieliniowe');