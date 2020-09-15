%Reconstruction characterisation V9
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
%V9: Determining resolutions and limits of reconstruction methods without
%    going throught the interferogram and phase extraction steps


close all

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
%% Creation of the 2D index profile to reconstruct
tic

index = zeros(N,N);

     %430
     
%Rectangular test object

%     width = 20;
%     height = 50;
%      
%      for i = N/2-(round(height/2)):N/2+(round(height/2))
%          for j = N/2-(round(width/2)):N/2+(round(width/2))
%              index(i,j) = n-1;
%          end
%      end 
%      



%Realistic test object


    %sigma = 25*(h-1)/(H-1)+40;
    %sigma = (2*(h-1))*70/(N-2) +10;
    sigma = 50
     
    
    
    for i = 1:N
        for j = 1:N
            dist = sqrt(((i-((N+1)/2))^2) + ((j-((N+1)/2))^2));
            %sigma = sqrt(d/2.5);
            %sigma = 7;
            if i > round((0.53*N))
                index(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(5)))*exp(-((((i-round(0.53*N))^2)/(2*10^2))^3));
                
            
            %if d < ((r/(L/(N-1))/N))*2*h + r/(L/(N-1))  %Circular index profile at each height
                %index(i,j) = ((1-n)/(((r/(L/(N-1))/N))*h + r/(L/(N-1))))*d + n;
            else
                %index(i,j) = n* exp(-(d/(2*sigma^2))^(-h*(8/N)+10) );
                
                index(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(5)));
                
                %layer(i,j) = (n-1)* exp(-(((dist^2)/(2*sigma^2))^(2)));
            %end
            
             end   
        end
    end

density_ref = (index*2*epsilon0)/(alpha*10^6);    %calculating density from index
    
figure
imagesc(density_ref)
colorbar
 xlabel('Z')
 ylabel('X')
 set(gca,'XTick',[], 'YTick', [])
title('Realistic reference image density profile')
% 
% figure
% imagesc(1+index(:,:,H/2))
% colorbar
% xlabel('Z')
% ylabel('X')
% set(gca,'XTick',[], 'YTick', [])
% %set(ax,'xticklabel',[])
% axis('xy')
% title('Index profile at a H/2 height')
% 
% figure
% imagesc(index(:,:,H))
% colorbar
% 
% figure
% imagesc(1+reshape(index(N/2,:,:),N,H)')
% colorbar
% xlabel('Z')
% ylabel('Y')
% set(gca,'XTick',[], 'YTick', [])
% axis('xy')
% title('Index profile at the middle of the jet (fixed x)')
% 
% figure
% imagesc(1+reshape(index(:,N/2,:),N,H)')
% colorbar
% xlabel('X')
% ylabel('Y')
% set(gca,'XTick',[], 'YTick', [])
% axis('xy')
% title('Index profile at the middle of the jet (fixed z)')



%figure
%imagesc(density(:,:,H/2))
%colorbar

%n = 1+ ((dens*10^6 * alpha ) / (2*epsilon0));
toc
%% Gas phase difference construction
tic

p = N;          %Number of rays for each angle
d = p-1;        %Distance between first and last ray
isDisp = 0; % No display
Xref = density_ref;
k = 25;             %Choose number of SIRT iterations
%sig_noise = 0.1 ;    %Choose sigma^2 of gaussian noise
sig_noise = 0:0.6:3;    %Choose range of sigmas

%angle_number = [2:8:66  150];  %Choose a range of numbers of used angles
%angle_number = [2 5 10];       %calculation with 2 angles then 5 then 10
angle_number =1;

%iterations_number = 1:20:500;  %k range

E = zeros(length(angle_number),3);  %Initializing error vector



%for ii = 1:length(angle_number)            %to plot error vs number of
                                            %angles used
                                            
for ii = 1:length(sig_noise)                %to plot error vs noise level
    
  %for ii = 1:length(iterations_number)     %to plot error vs nb of
                                            %iterations
    
ii

k = iterations_number(ii);

theta = 0:9:180;
%theta = 0:180/angle_number(ii):180;        %Angles for error vs number of
                                                                   %angles
                                                                   
%theta = [0 5 10 15 20 30 40 45 50 60 70 85 90 95 110 120 135 150 165 180];
%theta = [10];


[A] = paralleltomo(N,theta,p,d);

delta_phi_gas = zeros(N,length(theta));       %dimensions: H x N x Number of angles
    
image = index(:);
delta_phi_gas = reshape(A*image,N,length(theta));
%B = A*image;
%delta_phi_gas = b;

%imagesc(b)
%imagesc(delta_phi_gas)

delta_phi_gas = delta_phi_gas *(L/N)*(2*pi/lambda);

% figure
% imagesc(delta_phi_gas)
% %axis('xy')
% colorbar

% Reconstruction

B_gen = delta_phi_gas(:);          %Generated B

noise_b = imnoise(zeros(size(B_gen)),'gaussian',0,sig_noise);
B = B_gen + noise_b(:);

%Landweber Method
Xland = abs(landweber(A,B,k));
index_land = Xland./((L/N)*(2*pi/lambda));
density_land = (index_land*2*epsilon0)/(alpha*10^6);
    
%Cimmino Method
Xcim = abs(cimmino(A,B,k));
index_cim = Xcim./((L/N)*(2*pi/lambda));
density_cim = (index_cim*2*epsilon0)/(alpha*10^6);


%CAV Method
Xcav = abs(cav(A,B,k));
index_cav = Xcav./((L/N)*(2*pi/lambda));
density_cav = (index_cav*2*epsilon0)/(alpha*10^6);

E(ii,1) = norm(density_land(:)-density_ref(:));
E(ii,2) = norm(density_cim(:)-density_ref(:));
E(ii,3) = norm(density_cav(:)-density_ref(:));
       
end 
    

E = E/(norm(density_ref(:)))

figure
plot(sig_noise,E)
xlabel('k Number of iterations')
ylabel('Relative error')
legend('Landweber','Cimmino','CAV','location','best')
title('Relative error using 3 methods vs number of SIRT iterations')

toc
%%

    %Display Reference Image
figure
imagesc(reshape(Xref,N,N))
xlabel('Z')
ylabel('X')
set(gca,'XTick',[], 'YTick', [])
colorbar
title("Reference image density profile" )
    
    %Display Landweber reconstruction
        figure
        imagesc(reshape(density_land,N,N))
        axis image on
        colorbar
        title('Landweber')
    
    %Display Cimmino reconstruction
        figure
        imagesc(reshape(density_cim,N,N))
        axis image off
        colorbar
        title('Cimmino')
        %%
        
    %Display CAV reconstruction
        thresh = 0;          %Threshold to apply in % of max
        thresh_filter = zeros(N*N,1);
        thresh_filter(density_cav >thresh*max(density_cav)/100) = 1;
        
        figure
        imagesc(reshape(density_cav.*thresh_filter,N,N))
        %imagesc(reshape(density_cav,N,N))
        colorbar
        xlabel('Z')
        ylabel('X')
        set(gca,'XTick',[], 'YTick', [])
        title("Reconstruction using CAV and {\sigma^2_{noise}} ="+sig_noise+" Threshold = "+thresh + "%")

        E_filtered = norm(density_cav.*thresh_filter - density_ref(:))/norm(density_ref(:))
        
        