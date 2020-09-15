%3D gas cone interferometry V3
%15 06 2020 Start

%This script simulates a 2D Mach-Zehnder interferometer in which one path
%contains a disk jet of gas of index n(x,y,z). The goal is to plot the resulting
%intensity on a screen which corresponds to the output of the device.

%V1: Gas index profile is a definite function of space, analytical model
%V2: Gas index profile is an arbitrary discrete function of space
%V3: Introducing height, 2D interferogram
%



close all



lambda = 632.8 * 10^-9;      %HeNe wavelength in m
%n = 1.00027751;             %Optical index of Argon gas @632.8nm normal conditions
n = 1.00001751;
r = 10 * 10^-3;                      %Radius of the gas flow in m
I0 = 1;                     %Incident intensity
alpha_deg = 0.01;
alpha_rad = alpha_deg *pi/180;
L = 60*10^-3;   %Length of gas jet analysis zone
N = 1001;      %Corresponds to 1000 intervals 

y = linspace(-L/2,L/2,N);



for h = 1:N
    
    index = ones(N,N);

    for i = 1:N
        for j = 1:N
            if sqrt(((i-((N+1)/2))^2) + ((j-((N+1)/2))^2))< ((r/(L/(N-1))/N))*h + r/(L/(N-1))  %Circular index profile at each height
                index(i,j) = n;
            end
        end
    end


delta_phi_gas(N+1-h,:) = (L/N)*(2*pi/lambda) .* sum(index-ones(N,N),2)';

end

delta_phi_tilt = ((2*pi/lambda) * alpha_rad * ones(N,1) * y)'; %Phase shift due to tilt of angle alpha between the two wavefronts
i = lambda/alpha_rad;   %Fringe spacing of tilt fringes in m

I_screen = I0*(1+cos(delta_phi_gas+delta_phi_tilt));
 


%% Generation of interferograms

close all
% figure
% imshow(delta_phi_tilt)
% colorbar

figure
imshow(I_screen)
title('Interferogram image associated with a gas cone')
xlabel('x position on screen')
ylabel('y position on screen')
%%
imwrite(I_screen, 'Cone_10_20.jpg');







%%