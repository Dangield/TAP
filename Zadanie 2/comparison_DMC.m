clear;
close all;
clc;

addpath('./abstraction')
addpath('./classes')
addpath('./..')

consts
load("data/setPoints.mat")

load("data/numericalResult.mat")
un = u;
yn = y;

load("data/analyticalResult.mat")

ua = u;
ya = y;

titles = ["C_A", "T"];
figure;
for i = 1:2
    subplot(2, 1, i);
        hold on;
        stairs(setPoints(i, :), 'r');
		stairs(ya(i, :), 'b');
		stairs(yn(i, :), 'm');
		title(titles(i));
		legend("Wartosci zadane", "DMC Analityczny", "DMC numeryczny")
end


titles = ["C_Ain", "F_c"];
figure;
for i = 1:2
    subplot(2, 1, i);
        hold on;
		stairs(ua(i, :), 'b');
		stairs(un(i, :), 'm');
		title(titles(i));
		legend("DMC Analityczny", "DMC numeryczny")
end

ea = (ya - setPoints)*(ya - setPoints)';
ea = [ea(1, 1) ea(2, 2)]

en = (yn - setPoints)*(yn - setPoints)';
en = [en(1, 1) en(2, 2)]

ea - en