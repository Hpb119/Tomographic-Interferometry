%2D gas disk interferometry V1
%12 06 2020 Start

%This script simulates a 2D Mach-Zehnder interferometer in which one path
%contains a disk jet of gas of index n. The goal is to plot the resulting
%intensity on a screen which corresponds to the output of the device.

%V1: Gas index profile is a definite function of space, analytical model

%Working version 15 06 2020

close all



lambda = 632.8 * 10^-9;      %HeNe wavelength in m
n = 1.00027751;             %Optical index of Argon gas @632.8nm normal conditions
r = 10 * 10^-3;                      %Radius of the gas flow in m
I0 = 1;                     %Incident intensity
alpha_deg = 0.01;
alpha_rad = alpha_deg *pi/180;

x = linspace(-0.05,0.05,10001);    %Creation of the 1D study points in m

e_x = real(2*sqrt(r^2-x.^2));      %Thickness profile of gas jet

delta_phi_gas = (2*pi/lambda).*e_x*(n-1); %Phase shift induced by the presence of n index gas

delta_phi_tilt = (2*pi/lambda) * alpha_rad * x; %Phase shift due to tilt of angle alpha between the two wavefronts
i = lambda/alpha_rad;   %Fringe spacing of tilt fringes in m

I_screen = I0*(1+cos(delta_phi_gas+delta_phi_tilt)); %Resulting output intensity


figure
plot(x,e_x);     %Plotting the profile
movegui('northwest')
axis equal

figure
plot(x,delta_phi_tilt + delta_phi_gas);
movegui('north')

figure
plot(1000*x,I_screen);
xlabel('x position on screen in mm')
ylabel('Output intensity on screen')
title('1D interferogram generated by a centered 10mm radius gas jet disk and 0.01� tilt')
movegui('northeast')