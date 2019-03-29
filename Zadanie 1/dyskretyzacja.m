close all
clear
addpath("../")

consts
linearyzacja;
clc

c2d(transferFunction, 0.1)
figure
step(ans)

c2d(transferFunction, 0.5)
figure
step(ans)