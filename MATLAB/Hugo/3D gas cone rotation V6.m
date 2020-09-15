
%3D gas cone interferometry V3
%
%This script simulates a 2D Mach-Zehnder interferometer in which one path
%contains a disk jet of gas of index n(x,y,z). The goal is to plot the resulting
%intensity on a screen which corresponds to the output of the device.

%V1: Gas index profile is a definite function of space, analytical model
%V2: Gas index profile is an arbitrary discrete function of space
%V3: Introducing height, 2D interferogram
%V4: Rectangular rather than conic profile
%V5: Non uniform index
%V6: Rotation of angle theta + asymmetry
%



close all

tic

lambda = 632.8 * 10^-9;      %HeNe wavelength in m

dens = 5 * 10^18;   %Particle density in atoms/cm^3
alpha = 18.52 * 10^(-41);       %Polarizability of Argon in F m^2
epsilon0 = 8.85 * 10^-12;       % Permittivity of free space in F m^-1

n = 1+ ((dens*10^6 * alpha ) / (2*epsilon0));

%n = 1.00027751;             %Optical index of Argon gas @632.8nm normal conditions

theta = 0;

r = 10 * 10^-3;                      %Radius of the gas flow in m
I0 = 1;                     %Incident intensity
alpha_deg = 0.01;
alpha_rad = alpha_deg *pi/180;
L = 60*10^-3;   %Length of gas jet analysis zone
N = 1001;      %Corresponds to 1000 intervals 

y = linspace(-L/2,L/2,N);

delta_phi_gas = zeros(N,N);

for h = 1:N/2
    
    index = zeros(N,N);
    %sigma = 50;
    sigma = (2*(h-1))*70/(N-2) +70;
    
    for i = 430:N
        for j = 1:N
            d = sqrt(((i-((N+1)/2))^2) + ((j-((N+1)/2))^2));
            %sigma = sqrt(d/2.5);
            %sigma = 7;
            
            if d < ((r/(L/(N-1))/N))*2*h + r/(L/(N-1))  %Circular index profile at each height
                %index(i,j) = ((1-n)/(((r/(L/(N-1))/N))*h + r/(L/(N-1))))*d + n;
                %index(i,j) = n* exp(-(d/(2*sigma^2))^(-h*(8/N)+10) );
                index(i,j) = (n-1)* exp(-(((d^2)/(2*sigma^2))^(-2*(h-1)*(8/(N-2))+10)));
            end
        end
    end

index((N+1)/2,(N+1)/2) = n-1;

index = imrotate(index,theta,'crop');

delta_phi_gas(N+1-(2*h),:) = (L/N)*(2*pi/lambda) .* sum(index,2)';
delta_phi_gas(N+1-(2*h-1),:) = (L/N)*(2*pi/lambda) .* sum(index,2)';

end

delta_phi_tilt = ((2*pi/lambda) * alpha_rad * ones(N,1) * y)'; %Phase shift due to tilt of angle alpha between the two wavefronts
i = lambda/alpha_rad;   %Fringe spacing of tilt fringes in m

I_screen = I0*(1+cos(delta_phi_gas+delta_phi_tilt));

toc

%% Generation of interferograms

close all
% figure
% imshow(delta_phi_tilt)
% colorbar

figure
imshow(I_screen)

%imshow(I0*(1+cos(delta_phi_tilt))) %ref
title('Interferogram image associated with a gas cone')
xlabel('x position on screen')
ylabel('y position on screen')
%%
imwrite(I0*(1+cos(delta_phi_tilt)), 'Ref_SupergaussianV6.jpg');     %Ref

imwrite(I_screen, 'Cone_SupergaussianV6.jpg');





%%