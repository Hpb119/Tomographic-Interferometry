%Tomographic_interferometry_Final
%15/09/2020
%
%This script simulates a 2D Mach-Zehnder interferometer in which one path
%contains a disk jet of gas of index n(x,y,z). The goal is to plot the 
%resulting intensity on a screen which corresponds to the output of the 
%device. We then extract the corresponding phase difference and reconstruct
%a 3D image using multiple viewing angles.

%V1: Gas index profile is a definite function of space, analytical model
%V2: Gas index profile is an arbitrary discrete function of space
%V3: Introducing height, 2D interferogram
%V4: Rectangular rather than conic profile
%V5: Non uniform index
%V6: Rotation of angle theta + asymmetry
%V7: Sinogram from phase maps of different angles
%V8: Separating steps + better projection generation + functionnal
%reconstruction



%%%%%%%Execute sections one by one

close all
clear

%% Constant parameters
tic

lambda = 632.8 * 10^-9;      %HeNe wavelength in m

dens = 5 * 10^18;   % Max particle density in atoms/cm^3
alpha = 18.52 * 10^(-41);       % Polarizability of Argon in F m^2
epsilon0 = 8.85 * 10^-12;       % Permittivity of free space in F m^-1

n = 1+ ((dens*10^6 * alpha ) / (2*epsilon0)); % Index function of density
%n = 1.00027751;             %Optical index of Argon gas @632.8nm normal conditions

r = 10 * 10^-3;                 %Radius of the gas flow in m
I0 = 1;                         %Incident intensity
alpha_deg = 0.01;
alpha_rad = alpha_deg *pi/180;
L = 60*10^-3;   %Length of gas jet analysis zone in m
N = 500;        %Number of pixels along horizontal axes
H = 200;        %Number of pixels in height
l_ratio = L/N;  %Length corresponding to 1 pixel
h_ratio = L/H;  %Length corresponding to 1 pixel in height

toc
%% Creation of the 3D index profile
%Longest simulation step - ~12 seconds for H = 200, N = 500

tic

index = zeros(N,N,H);   %Initialisation of the 3D index profile

for h = 1:H             %Attributing values height by height
    
    h
    sigma = 25*(h-1)/(H-1)+40;      %Supergaussian varying width parameter
    %sigma = (2*(h-1))*70/(N-2) +10;
     %430
    layer = zeros(N,N);
    for i = 1:N
        for j = 1:N
            dist = sqrt(((i-((N+1)/2))^2) + ((j-((N+1)/2))^2)); %distance from center of matrix in pixels
            
            if i > round((0.53*N))  %Creating the asymmetry
                layer(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(-(h-1)*(8/(H-1))+10)))*exp(-((((i-round(0.53*N))^2)/(2*10^2))^3));
                
            
            %if d < ((r/(L/(N-1))/N))*2*h + r/(L/(N-1))  %Circular index profile at each height
                %index(i,j) = ((1-n)/(((r/(L/(N-1))/N))*h + r/(L/(N-1))))*d + n;
            else
                %index(i,j) = n* exp(-(d/(2*sigma^2))^(-h*(8/N)+10) );
                
                %Circular supergaussian profile at each height
                layer(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(-(h-1)*(8/(H-1))+10)));
                
                %layer(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(2)));
            %end
            
             end   
        end
    end
    index(:,:,h) = layer(:,:);
end

%Displaying profile views

figure
imagesc(index(:,:,1))
colorbar
title('Index profile at the nozzle')

figure
imagesc(1+index(:,:,H/2))
colorbar
xlabel('Z')
ylabel('X')
set(gca,'XTick',[], 'YTick', [])
%set(ax,'xticklabel',[])
axis('xy')
title('Index profile at a H/2 height')

figure
imagesc(index(:,:,H))
colorbar

figure
imagesc(1+reshape(index(N/2,:,:),N,H)')
colorbar
xlabel('Z')
ylabel('Y')
set(gca,'XTick',[], 'YTick', [])
axis('xy')
title('Index profile at the middle of the jet (fixed x)')

figure
imagesc(1+reshape(index(:,N/2,:),N,H)')
colorbar
xlabel('X')
ylabel('Y')
set(gca,'XTick',[], 'YTick', [])
axis('xy')
title('Index profile at the middle of the jet (fixed z)')


density = (index*2*epsilon0)/(alpha*10^6);
%figure
%imagesc(density(:,:,H/2))
%colorbar

%n = 1+ ((dens*10^6 * alpha ) / (2*epsilon0));
toc
%% Gas phase difference construction
tic
theta = 0:10:180;   %Viewing angles to consider
p = N;          %Number of rays for each angle
d = p-1;        %Distance between first and last ray
isDisp = 0; % No display

[A] = paralleltomo(N,theta,p,d);            %Use of the AIR II function to generate matrix A
%A = get_or_apply_system_matrix(N,theta,p,d,isDisp);
delta_phi_gas = zeros(H,N,length(theta));       %dimensions: H x N x Number of angles


for h = 1:H                               %For each height
    image = index(:,:,h);                   %Generation of the cross section at the height
    image = image(:);                       % Building a column vector
    b = reshape(A*image,N,length(theta));   %Calculating the projections for all angles for this cross section
    %B = A*image;
    delta_phi_gas(h,:,:) = b;               %Storing projections
end
%imagesc(b)
%imagesc(delta_phi_gas)

delta_phi_gas = delta_phi_gas *(L/N)*(2*pi/lambda); %Converting into phase

disp_angle_nb = 6;                                  %index of the angle in 'theta' which you want to display in figures

max_phase = max(max(delta_phi_gas(:,:,disp_angle_nb))); %Calculations to scale the colorbar
phase_range = 0:pi:fix(max_phase/pi)*pi;

figure
imagesc(delta_phi_gas(:,:,disp_angle_nb))       %Display of a phase map at a given angle
axis('xy')
colorbar
colorbar('Ticks', phase_range, 'TickLabels',['0   ';' \pi';'2\pi';'3\pi';'4\pi';'5\pi';'6\pi';'7\pi'])
xlabel('X')
ylabel('Y')
set(gca,'XTick',[], 'YTick', [])
title("Phase shift map for a " + theta(disp_angle_nb) + " degrees viewing angle")


disp_height = 0.9*H;    %Display height

max_phase2 = max(max(delta_phi_gas(disp_height,:,:)));
phase_range2 = 0:pi:fix(max_phase2/pi)*pi;

figure
imagesc(reshape(delta_phi_gas(disp_height,:,:),N,length(theta)))

%set(gca,'XTick',1:10:61,'XTickLabel',0:30:180,'YTick',[])  % 
set(gca,'XTick',1:1:length(theta),'XTickLabel',theta,'YTick',[])

colorbar
colorbar('Ticks', phase_range2, 'TickLabels',['0   ';' \pi';'2\pi';'3\pi';'4\pi';'5\pi';'6\pi';'7\pi'])
xlabel('Viewing angle (deg)')
ylabel('X')
title("Generated sinogram using " + length(theta) + " angles at height H/2") 

toc


%% Interferogram construction
tic

disp_angle_nb = 5;   %Number of the angle to display the interferogram of
y = linspace(-L/2,L/2,H);

noisy = true;       %add Gaussian white noise

noise_m = 0;    %Mean value of noise
noise_v = 0.1;  %Variance of noise

noise = zeros(H,N,length(theta));


delta_phi_tilt_2D = ((2*pi/lambda) * alpha_rad * ones(N,1) * y)'; %Phase shift due to tilt of angle alpha between the two wavefronts
delta_phi_tilt_3D = delta_phi_tilt_2D;
for th = 1:length(theta)-1
    
    delta_phi_tilt_3D = cat(3,delta_phi_tilt_3D,delta_phi_tilt_2D);
end

fringe_spacing = lambda/alpha_rad;   %Fringe spacing of tilt fringes in m

I_screen_simu = I0*(1+cos(delta_phi_gas+delta_phi_tilt_3D));    %Perfect intensity on the screen

if noisy
    noise = imnoise(zeros(H,N,length(theta)),'gaussian',noise_m,noise_v);   %add gaussian noise of mean noise_m and sigma^2 = noise_v 
end

I_screen = I_screen_simu + noise;   %Intensity on the screen + noise


figure
imagesc(I_screen(:,:,disp_angle_nb))
xlabel('X')
ylabel('Y')
set(gca,'XTick',[], 'YTick', [])
colorbar
h=colorbar;
t=get(h,'Limits');
colorbar('Ticks', [t(1) 1 t(2)], 'TickLabels',[0 1 t(2)])
axis('xy')
colormap gray
title("Interferogram for a " + theta(disp_angle_nb) + "-degree viewing angle. \sigma_{noise}^2 = "+noise_v)


toc
%% Phase extraction
%Adapted from https://github.com/jasmcole/interpret
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure
% imshow(I_screen)
% 
% title("Interferogram for a " + theta + " degree angle")
% xlabel('x position on screen')
% ylabel('y position on screen')
tic
extracted_phase = zeros(H,N,length(theta));

pick_fourier_region = false;     %picking new fourier region

%Default Fourier region            coordinates of a corner + width + height
x2 = 65.5;
y2 = 139.1;
w = 27.95;
v = 222.7;


for ang = 1:length(theta)
%ang = 13;

data = double(imrotate(I_screen(:,:,ang),90));
intref = double(imrotate(noise(:,:,ang) + I0*(1+cos(delta_phi_tilt_3D(:,:,ang))),90));

fftim = fft2(data);
fftref = fft2(intref);

%     
if ang == 1
    if pick_fourier_region
        figure
        imagesc(log(abs(fftshift(fftim))))
        
        set(gca,'XTick',[], 'YTick', [])
        title('Visualisation of the 2D FFT of an interferogram')
        zone = getrect
        x2 = zone(1);
        y2 = zone(2);
        w = zone(3);
        v = zone(4);
        
        close 
    end  
end



fftim = fftshift(fftim);
fftref = fftshift(fftref);

pow = 16;
[Ny Nx] = size(fftim);
[xg yg] = meshgrid(1:Nx, 1:Ny);
mask = exp( -((xg - x2 - w/2).^pow)/((w/2)^pow)).*exp( -((yg - y2 - v/2).^pow)/((v/2)^pow));

newfftim = mask.*fftim;
newfftref = mask.*fftref;


newfftim = double(newfftim);
newfftref = double(newfftref);
newfftim = fftshift(newfftim);
newfftref = fftshift(newfftref);

datafiltered = ifft2(newfftim);
reffiltered = ifft2(newfftref);

if(sum(sum(isnan(reffiltered))) > 0)
    reffiltered = ones(size(datafiltered));
end
if (max(max(intref)) == min(min(intref)))
    reffiltered = ones(size(datafiltered));
end
phase = imrotate(angle(datafiltered./reffiltered),-90);
mag = abs(datafiltered./reffiltered);

%Unwrapping
phase = flipud(unwrap(flipud(phase)));
phase = imrotate(phase, 90);
phase = unwrap(phase);
phase = imrotate(phase, -90);

extracted_phase(:,:,ang) = abs(phase);       %H x N x length(theta) matrix 
                                        %containing extracted phase info
end
% figure
 %imagesc(phase); colorbar
% title("Phase map for a " + theta + " degree angle")

max_phase2 = max(max(extracted_phase(disp_height,:,:)));
phase_range3 = 0:pi:fix(max_phase2/pi)*pi;

figure
imagesc(reshape(extracted_phase(disp_height,:,:),N,length(theta)))

set(gca,'XTick',1:1:length(theta),'XTickLabel',theta,'YTick',[])

colorbar
colorbar('Ticks', phase_range3, 'TickLabels',['0   ';' \pi';'2\pi';'3\pi';'4\pi';'5\pi';'6\pi';'7\pi'])
xlabel('Viewing angle (deg)')
ylabel('X')
title('Extracted sinogram at height 0.9H')




max_phase = max(max(extracted_phase(:,:,disp_angle_nb)));
phase_range = 0:pi:fix(max_phase/pi)*pi;

figure
imagesc(extracted_phase(:,:,disp_angle_nb))
set(gca,'XTick',[], 'YTick', [])
axis('xy')
colorbar
colorbar('Ticks', phase_range, 'TickLabels',['0   ';' \pi';'2\pi';'3\pi';'4\pi';'5\pi';'6\pi';'7\pi'])
xlabel('X')
ylabel('Y')
title("Retrieved phase for a " + theta(disp_angle_nb) + "-degree viewing angle")

%Reconstructed vs real phase difference

% figure
% retrieval_error = norm(extracted_phase(:,:,disp_angle_nb) - delta_phi_gas(:,:,disp_angle_nb),'fro')
% mean(mean(abs(extracted_phase(:,:,disp_angle_nb) - delta_phi_gas(:,:,disp_angle_nb))))
% %imagesc(retrieval_error)
% axis('xy')
% colorbar



toc

%% Reconstruction
%Using function and structures from thesis reference AIR II.
%%%%%%%%%%%%
tic


height = round(120);

B = reshape(extracted_phase(height,:,:),N,length(theta));
B = B(:);

nMethod = 3;
k = 100;

Xref = density(:,:,height);
Xsart = abs(sart(A,B,k));
index_sart = Xsart./((L/N)*(2*pi/lambda));
density_sart = (index_sart*2*epsilon0)/(alpha*10^6);

figure
subplot(double(nMethod>1)+1,2,1)
imagesc(reshape(Xref,N,N))
axis image off
colorbar
title("Reference at height "+ height )

subplot(double(nMethod>1)+1,2,2)
imagesc(reshape(density_sart,N,N))
axis image off
colorbar
title("SART" )


if nMethod>1
    
Xdrop = abs(drop(A,B,k));
index_drop = Xdrop./((L/N)*(2*pi/lambda));
density_drop = (index_drop*2*epsilon0)/(alpha*10^6);
 
 % Show the DROP solution.
 subplot(double(nMethod>1)+1,2,3)
 imagesc(reshape(density_drop,N,N))
 axis image off
 colorbar
 title('DROP')
    if nMethod>2
Xland = abs(landweber(A,B,k));
index_land = Xland./((L/N)*(2*pi/lambda));
density_land = (index_land*2*epsilon0)/(alpha*10^6);
 
 % Show the DROP solution.
 subplot(double(nMethod>1)+1,2,4)
 imagesc(reshape(density_land,N,N))
 axis image off
 colorbar
 title('Landweber')
    end
end 
sgtitle('Reference vs reconstructed density profiles')

toc


%%
%Saving images

imwrite(I0*(1+cos(delta_phi_tilt)), "Ref_"+theta+"deg_V7.jpg");     %Ref

imwrite(I_screen, "Cone_"+theta+"deg_V7.jpg"); 


