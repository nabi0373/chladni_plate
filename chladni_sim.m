%% Chladni Plate Analytic Solution Simulator 
%clear all;
%close all;

compute_circle_model = false; %only set to true if circle data isn't already cached
numcirc_to_find = 6; %number of discrete frequencies/wavelengths to find

compute_annulus_model = false; %only set to true if annulus data isn't already cached
numann_to_find = 219; %number of discrete frequencies/wavelengths to find

compute_driver_decomposition = false;
save_driver_decomposition = false;
load_driver_decomposition = true;
csv_filename = 'flat_driving.csv';

run_chladni_animation_2d = true;
run_chladni_animation_3d = false;

%% Definition of variables
Ri = 0.0021; %in meters
Ro = 0.12; % in meters
h = 0.00095; %in meters
rho = 2700; %in kg/m^3
PR = 0.33; %Poisson Ratio
YM = 68.9e9; %Young's Modulus in kg/(m*s^2)
D = YM*h^3/(12*(1-PR^2));
beta = rho*h/D;

%% Definition of Bessel Functions
syms lambda r t
J0o = besselj(0,lambda*Ro);
J0i = besselj(0,lambda*Ri);
J1o = besselj(1,lambda*Ro);
J1i = besselj(1,lambda*Ri);
J2o = besselj(2,lambda*Ro);
J2i = besselj(2,lambda*Ri);
J3o = besselj(3,lambda*Ro);
J3i = besselj(3,lambda*Ri);

Y0o = bessely(0,lambda*Ro);
Y0i = bessely(0,lambda*Ri);
Y1o = bessely(1,lambda*Ro);
Y1i = bessely(1,lambda*Ri);
Y2o = bessely(2,lambda*Ro);
Y2i = bessely(2,lambda*Ri);
Y3o = bessely(3,lambda*Ro);
Y3i = bessely(3,lambda*Ri);

I0o = besseli(0,lambda*Ro);
I0i = besseli(0,lambda*Ri);
I1o = besseli(1,lambda*Ro);
I1i = besseli(1,lambda*Ri);
I2o = besseli(2,lambda*Ro);
I2i = besseli(2,lambda*Ri);
I3o = besseli(3,lambda*Ro);
I3i = besseli(3,lambda*Ri);

K0o = besselk(0,lambda*Ro);
K0i = besselk(0,lambda*Ri);
K1o = besselk(1,lambda*Ro);
K1i = besselk(1,lambda*Ri);
K2o = besselk(2,lambda*Ro);
K2i = besselk(2,lambda*Ri);
K3o = besselk(3,lambda*Ro);
K3i = besselk(3,lambda*Ri);

%% Circle Model
if compute_circle_model
    circle_matrix = [(1/2)*(-J0o+J2o),(1/2)*(I0o+I2o); ...
                     (1/4)*(3*J1o-J3o),(1/4)*(3*I1o+I3o)];
    det_circle = det(circle_matrix);
    f = matlabFunction(det_circle);
    circle_lambdas = [];
    tolerance = 1e-10;

    i = 0; %Increasing variable for zero guesses
    while length(circle_lambdas) < numcirc_to_find + 1
        i = i+1;
        zero = fzero(@(lambda) f(lambda),i);
    
        if ~checkInZeros(circle_lambdas,zero,tolerance)
            circle_lambdas = [circle_lambdas,zero];
        end
    end

    circle_lambdas = circle_lambdas(abs(circle_lambdas) > tolerance); %Remove root at lambda = 0 if it got added
    circle_omegas = circle_lambdas.^2 ./sqrt(beta);
    circle_frequencies = circle_omegas/(2*pi);

    % Next, find the nullspace of circle_matrix for each lambda in
    % circle_lambdas and find the symbolic equation for each solution
    circle_nullspaces = [];
    circle_equations = [];

    for i = circle_lambdas
        nullspace = null(feval(matlabFunction(circle_matrix),i))';
        if nullspace(1) < 0
            nullspace = nullspace * -1;
        end
        circle_nullspaces = [circle_nullspaces; nullspace];
        circle_equations = [circle_equations,nullspace(1)*besselj(0,i*r)+nullspace(2)*besseli(0,i*r)];
    end
end

%% Annulus Model
if compute_annulus_model
    annulus_matrix = [J0i,Y0i,I0i,K0i; ...
                      -J1i,-Y1i,I1i,-K1i; ...
                      (1/2)*(-J0o+J2o),(1/2)*(-Y0o+Y2o),(1/2)*(I0o+I2o),(1/2)*(K0o+K2o); ...
                      (1/4)*(3*J1o-J3o),(1/4)*(3*Y1o-Y3o),(1/4)*(3*I1o+I3o),(-1/4)*(3*K1o+K3o)];

    det_annulus = det(annulus_matrix);
    f = matlabFunction(det_annulus);
    annulus_lambdas = [];
    tolerance = 1e-10;

    i = 10; %Increasing variable for zero guesses
    while length(annulus_lambdas) < numann_to_find
        i = i+1;
        disp(length(annulus_lambdas));
        disp(i);
        disp(f(i));
        zero = fzero(@(lambda) f(lambda),i);
    
        if ~checkInZeros(annulus_lambdas,zero,tolerance)
            annulus_lambdas = [annulus_lambdas,zero];
        end
    end

    annulus_lambdas = annulus_lambdas(abs(annulus_lambdas) > tolerance); %Remove root at lambda = 0 if it got added
    annulus_omegas = annulus_lambdas.^2 ./sqrt(beta);
    annulus_frequencies = annulus_omegas/(2*pi);

    % Next, find the nullspace of circle_matrix for each lambda in
    % circle_lambdas and find the symbolic equation for each solution
    annulus_nullspaces = [];
    annulus_equations = [];

    for i = annulus_lambdas
        nullspace = null(feval(matlabFunction(annulus_matrix),i))';
        if nullspace(1) > 0
            nullspace = nullspace * -1;
        end
        annulus_nullspaces = [annulus_nullspaces; nullspace];
        annulus_equations = [annulus_equations,nullspace(1)*besselj(0,i*r)+nullspace(2)*bessely(0,i*r)+nullspace(3)*besseli(0,i*r)+nullspace(4)*besselk(0,i*r)];
    end
end

%% Driven Boundary Condition Model
if compute_driver_decomposition
    epsilon = (Ro-Ri)/100000;
    Rend = Ri+epsilon;
    A = 1; %1mm motor displacement
    driving_func = piecewise(r < Ri+epsilon,A,r >= Ri+epsilon,0); % F(r,theta) -- can be changed to other profiles
                                                                  % This is assuming periodic forcing 
    
    f = 0.991*annulus_frequencies(4); %Set sinusoidal driving frequency here
    ohm = 2*pi*f;
    Tcycle = 1/f;                                                          
                                                                  
    q = zeros(1,numann_to_find); %initialize empty array for coefficients of bessel basis functions
    driving_func_approx = 0*r; %approximate partial sum of driving function from bessel basis

    for i = 1:numann_to_find
    
        q(i) = double(int(driving_func*annulus_equations(i)*r,r,Ri,Ro)/int((annulus_equations(i)^2)*r,r,Ri,Ro));
        
        if save_driver_decomposition
            dlmwrite(csv_filename,q(i),'delimiter',',','-append');
        end
    
        driving_func_approx = driving_func_approx + annulus_equations(i)*q(i);
        
        %Print out the progress
        if mod(i,floor(numann_to_find/100)) == 0
            disp(sprintf('%3.2f%% done!',100*i/numann_to_find));
        end
    end
end

if load_driver_decomposition
    q = csvread('epsilon=0_00001.csv');
end

usolution = 0*t;
for i = 1:numann_to_find
    if f/annulus_frequencies(i) > 0.99 && f/annulus_frequencies(i) < 1/0.99 %if driving frequency near a resonant frequency
        usolution = usolution - (q(i)/(2*ohm))*t*cos(ohm*t)*annulus_equations(i); 
    else
        usolution = usolution + (q(i)/((1/beta)*annulus_lambdas(i)^4 - ohm^2))*(sin(ohm*t)-(ohm*sqrt(beta)/annulus_lambdas(i)^2)*sin(((1/sqrt(beta))*annulus_lambdas(i)^2)*t))*annulus_equations(i);
    end
end

disp('Done Computing Solution!');

%% Chladni Animation - Following https://www.mathworks.com/help/matlab/creating_plots/animating-a-surface.html

if run_chladni_animation
    u = matlabFunction(usolution); %spatial function to be animated

    %Define the polar grid
    rad_steps = 100;
    rarr = Ri:(Ro-Ri)/rad_steps:Ro;
    theta_steps = 50;
    theta = 0:2*pi/theta_steps:2*pi;
    [rad,theta] = meshgrid(rarr,theta); 

    %Convert to cartesian coordinates
    x = rad.*cos(theta);
    y = rad.*sin(theta);

    scaling_factor = 5e14;
    z = scaling_factor*u(rad,0);

    c = 128*ones(theta_steps,rad_steps); %RGB grey is (128,128,128)

    data = data*10;

    if run_chladni_animation_3d
        figure(1);
        s = surf(x,y,z,c);
        light               % add a light
        lighting gouraud    % preferred lighting for a curved surface
        axis equal off     % set axis equal and remove axis
        view(45,30)         % set viewpoint
        camzoom(1)        % zoom into scene
    end
    
    
    num_cycles = 30;
    num_steps_per_cycle = 100;
    tsteps = num_cycles*num_steps_per_cycle;
    tstart = 0;
    tend = tstart + Tcycle*num_cycles;
    k = (tend-tstart)/tsteps;
    tarr = tstart:k:tend;
    
    tdelay = 0.05;
    
    %pre-compute data at each time
%     data = zeros(51,101,length(tarr));
%     func = matlabFunction(usolution);
%     for i = 1:length(tarr)
%         disp(i);
%         for j = 1:length(rarr)
%             %disp(j);
%             %data(:,j,i) = double(subs(usolution,[r,t],[rarr(j),tarr(i)]));
%             data(:,j,i) = scaling_factor*feval(func,rarr(j),tarr(i));
%         end
%     end
    
    %data = data/100;
    
    if run_chladni_animation_2d
        figure(2);
        p = plot(1:1:101,data(1,:,1));
    end
    
    for i = 1:length(tarr)
        disp(i);
        
        if run_chladni_animation_3d
            figure(1);
            znew = data(:,:,i);
            %axis([-1.2*Ro 1.2*Ro -1.2*Ro 1.2*Ro -30 30]);
            s.ZData = znew;
        end
        
        pause(tdelay);
        
        if run_chladni_animation_2d
            figure(2);
            unew = data(1,:,i);
            p.YData = unew;
            axis([0,120,-30,30]);
        end
        
    end
end

function inZeros = checkInZeros(zeros,zero,tolerance)
    for k = zeros
        if abs(k-zero) < tolerance || abs(zero) < tolerance || isnan(zero)
            inZeros = true;
            return;
        end
    end
    inZeros = false; 
end

