%% Chladni Plate Finite Difference Method Solver
%Author: Nathan Biesterfeld
%Lives at https://github.com/nabi0373/chladni_plate

%% Script Parameters (Change these)
compute_chladni_solution = false; % set to true if you want to run the FDM solver
                                  % FDM solver parameters are marked
                                  % with %%%%%%%%%%FDM_PARAM%%%%%%%%%%
save_to_csv = false; % save results of FDM solver in a csv file specified by csv_filename.
                     % FDM solver data only saved if compute_chladni_solution is also true
csv_filename = '1941Hz_rsteps=100_percycle=7000_cycles=100.csv';

load_from_csv = false; % loads the solution matrix in csv_filename into sol_matrix
                       % Either load_solution_from_csv must be true
                       % OR compute_chladni_solution must be true
                       % (new chladni FDM solution is computed and stored in sol_matrix)
                       % OR sol_matrix must be cached from previous
                       % running of the script.

                                 
run_chladni_animation_3d = false; %runs the animation with the data in sol_matrix
run_chladni_animation_2d = true;
animation_pause = 0.00005; % time (in seconds) between frames of animation

plot_maximum_amplitudes = true; %plots maximum amplitudes of u(r,t) vs. time. 
                                 % Will see wave packets if driving
                                 % frequency isn't one of the natural
                                 % frequencies of the plate
                                 
plot_with_analytic_solution = false;
bessel_mode = 4; %The mode of analytic solution to plot numeric solution with
num_cycles_driver = 50; % The number of cycles of driver to plot the numeric solution after
                       % Must be less than num_cycles (the number of cycles
                       % in sol_matrix)
                                 
%% Definition of Plate Variables - Specific to the circular plate in the ITLL
Ri = 0.0021; %in meters
Ro = 0.12; % in meters
h = 0.00095; %in meters
rho = 2700; %in kg/m^3
PR = 0.33; %Poisson Ratio
YM = 68.9e9; %Young's Modulus in kg/(m*s^2)
D = YM*h^3/(12*(1-PR^2));
beta = rho*h/D;

%% Main call to computeChladniSolution

%{ Parameter values for various resonances
    %Temporal Mode     rsteps     num_steps_per_cycle     driving frequency
    %     2              20               300                    318
    %     2              50               2000                   350
    %     2              100              7000                   362
    %     3              20               300                    863
    %     3              50               2000                   956
    %     3              100              7000                   991
    %     4              20               300                    1675
    %     4              50               2000                   1871
    %     4              100              7000                   1941
    %     5              20               300                    2736
    %     5              50               2000                   3087
    %     5              100              7000                   3208
    %     6              20               300                    4027
    %     6              50               2000                   4603
    %     6              100              7000                   4792
%}

rsteps = 100; %%%%%%%%%%FDM_PARAM%%%%%%%%%%
h = (Ro-Ri)/rsteps;
r = Ri:(Ro-Ri)/rsteps:Ro; 

f = 1941; %driving frequency %%%%%%%%%%FDM_PARAM%%%%%%%%%%
Tcycle = 1/f;
ohm = 2*pi*f; %angular driving frequency

num_cycles = 100; %%%%%%%%%%FDM_PARAM%%%%%%%%%%
num_steps_per_cycle = 7000; %%%%%%%%%%FDM_PARAM%%%%%%%%%%
tsteps = num_cycles*num_steps_per_cycle;
tstart = 0;
tend = Tcycle*num_cycles;

k = (tend-tstart)/tsteps;
t = tstart:k:tend;

uinit = zeros(1,length(r)); %%%%%%%%%%FDM_PARAM%%%%%%%%%%
uinitinit = zeros(1,length(r)); %%%%%%%%%%FDM_PARAM%%%%%%%%%%

high_precision = 0; %%%%%%%%%%FDM_PARAM%%%%%%%%%%
                    % Set equal to zero for low precision (faster) and equal to one for high precision (slower) 
                    % Should be set to zero unless you have some powerful hardware to crank up num_steps_per_cycle 
                    % and force stability               

%Driven Boundary Conditions (These can be set to anything if you want to try other boundary conditions)
start_bc0_bool = true;
start_bc0 = sin(ohm.*t);
end_bc0_bool = false;
end_bc0 = zeros(1,size(t,2));
start_bc1_bool = true;
start_bc1 = zeros(1,size(t,2));
end_bc1_bool = false;
end_bc1 = zeros(1,size(t,2));
start_bc2_bool = false;
start_bc2 = zeros(1,size(t,2));
end_bc2_bool = true;
end_bc2 = zeros(1,size(t,2));
start_bc3_bool = false;
start_bc3 = zeros(1,size(t,2));
end_bc3_bool = true;
end_bc3 = zeros(1,size(t,2));
start_bc4_bool = false;
start_bc4 = zeros(1,size(t,2));
end_bc4_bool = false;
end_bc4 = zeros(1,size(t,2));

if compute_chladni_solution
    sol_matrix = computeChladniSolution(r,h,t,k,high_precision,uinit,uinitinit,beta,start_bc0_bool,start_bc0,start_bc1_bool,start_bc1,start_bc2_bool,start_bc2,start_bc3_bool,start_bc3,start_bc4_bool,start_bc4, ...
                                        end_bc0_bool,end_bc0,end_bc1_bool,end_bc1,end_bc2_bool,end_bc2,end_bc3_bool,end_bc3,end_bc4_bool,end_bc4,save_to_csv,csv_filename,load_from_csv);
end

%% Other tasks defined by script parameters
if load_from_csv
    sol_matrix = csvread(csv_filename);
    sprintf('Loaded solution data into sol_matrix from %s', csv_filename);
end

if plot_maximum_amplitudes
    plotMaxAmplitude(sol_matrix,t);
end

if plot_with_analytic_solution
    time = Tcycle*num_cycles_driver;
    plotNumericAnalytic(sol_matrix,r,bessel_mode,time,t,Ri,Ro)
end

%% Solution Animation
if run_chladni_animation_3d
    theta_steps = 100;
    theta = 0:2*pi/theta_steps:2*pi;
    [rad,theta] = meshgrid(r,theta); 
    
    %Convert to cartesian coordinates
    x = rad.*cos(theta);
    y = rad.*sin(theta);
    
    scaling_factor = 0.01;
    

    c = 128*ones(theta_steps+1,rsteps+1); %RGB grey is (128,128,128)
    
    z = zeros(theta_steps+1,rsteps+1);
    for i = 1:1:theta_steps+1
        z(i,:) = scaling_factor*sol_matrix(1,:);
    end
    
    fig = figure(3);
    s = surf(x,y,z,c);
    set(gcf,'Position',[0 0 1200 1200])
    light               % add a light
    lighting gouraud    % preferred lighting for a curved surface
    axis equal off      % set axis equal and remove axis
    view(45,30)         % set viewpoint
    camzoom(1.9)        % zoom into scene
    
    j = 1;
    for i = 2:length(t)
        znew = zeros(theta_steps+1,rsteps+1);
        for j = 1:1:theta_steps+1
            znew(j,:) = scaling_factor*sol_matrix(i,:);
        end
        s.ZData = znew;
        
        pause(animation_pause);
    end  
end

if run_chladni_animation_2d
    figure(4);
    p = plot(r,sol_matrix(1,:));
    ab_sol_matrix = abs(sol_matrix);
    axis([Ri, Ro, -max(ab_sol_matrix(:)), max(ab_sol_matrix(:))]);
    xlabel('Radial Distance from Center of Plate (m)');
    ylabel('Transverse Displacement');
    
    for i = 1:length(t)
        p.YData = sol_matrix(i,:);
        pause(animation_pause);
    end
end

%% Function for plotting maximum amplitude of u(r,t) as a function of time t (wave packets when not at resonance)
function plotMaxAmplitude(sol_matrix,t)
    assert(size(sol_matrix,1) == size(t,2));
    maximums = zeros(1,size(t,2)); %allocate memory for array of maximum values
    
    for i = 1:size(t,2)
    maximums(i) = max(abs(sol_matrix(i,:)));
    end
    
    figure(1);
    hold on;
    plot(t,maximums);
    xlabel('Time t (seconds)');
    ylabel('Maximum Amplitude of Plate at time t');
    title('Maximum Amplitude of plate vs. Time');
    hold off;
end

%% Function for plotting u(r,t) as a function of r at time t with annulus model Bessel function solution
function plotNumericAnalytic(sol_matrix,rarr,bessel_mode,time,tarr,Ri,Ro)
    
    syms r
    annulus_equations = [(2499828031264309*besseli(0, (7576134126845691*r)/562949953421312))/9007199254740992 - (5142445241127811*besselj(0, (7576134126845691*r)/562949953421312))/18014398509481984 + (4449207923115027*besselk(0, (7576134126845691*r)/562949953421312))/9007199254740992 + (6962317237895933*bessely(0, (7576134126845691*r)/562949953421312))/9007199254740992, (4902546510282909*besselk(0, (1311228229413527*r)/35184372088832))/9007199254740992 - (5282738823935647*besselj(0, (1311228229413527*r)/35184372088832))/288230376151711744 - (4487530694001037*besseli(0, (1311228229413527*r)/35184372088832))/288230376151711744 + (7552996552939613*bessely(0, (1311228229413527*r)/35184372088832))/9007199254740992, (4809513442131485*besseli(0, (2291821827263821*r)/35184372088832))/9223372036854775808 - (2651254712484407*besselj(0, (2291821827263821*r)/35184372088832))/36028797018963968 + (2498140102323221*besselk(0, (2291821827263821*r)/35184372088832))/4503599627370496 + (3732541978984283*bessely(0, (2291821827263821*r)/35184372088832))/4503599627370496, (2554765239672393*besselk(0, (1619546541027817*r)/17592186044416))/4503599627370496 - (8178554445780131*besselj(0, (1619546541027817*r)/17592186044416))/72057594037927936 - (91839264879511*besseli(0, (1619546541027817*r)/17592186044416))/4611686018427387904 + (7346917714268039*bessely(0, (1619546541027817*r)/17592186044416))/9007199254740992, (7353289709208021*besseli(0, (8363597359947883*r)/70368744177664))/9444732965739290427392 - (2782840742965823*besselj(0, (8363597359947883*r)/70368744177664))/18014398509481984 + (5236782042362697*besselk(0, (8363597359947883*r)/70368744177664))/9007199254740992 + (28105923870633*bessely(0, (8363597359947883*r)/70368744177664))/35184372088832, (5374419030160887*besselk(0, (2561272375484567*r)/17592186044416))/9007199254740992 - (7018781093492115*besselj(0, (2561272375484567*r)/17592186044416))/36028797018963968 - (4645309340715517*besseli(0, (2561272375484567*r)/17592186044416))/151115727451828646838272 + (7011868723846265*bessely(0, (2561272375484567*r)/17592186044416))/9007199254740992, (736862584146197*besseli(0, (757794224040501*r)/4398046511104))/604462909807314587353088 - (4211425816834423*besselj(0, (757794224040501*r)/4398046511104))/18014398509481984 + (2759859147621001*besselk(0, (757794224040501*r)/4398046511104))/4503599627370496 + (3399570616383139*bessely(0, (757794224040501*r)/4398046511104))/4503599627370496, (354410266761915*besselk(0, (1750410356251233*r)/8796093022208))/562949953421312 - (4879194171532747*besselj(0, (1750410356251233*r)/8796093022208))/18014398509481984 - (3747595563240179*besseli(0, (1750410356251233*r)/8796093022208))/77371252455336267181195264 + (6559169562836935*bessely(0, (1750410356251233*r)/8796093022208))/9007199254740992, (17347*besseli(0, (7940607665237127*r)/35184372088832))/9007199254740992 - (2752628078196721*besselj(0, (7940607665237127*r)/35184372088832))/9007199254740992 + (5825216868567003*besselk(0, (7940607665237127*r)/35184372088832))/9007199254740992 + (1573601075397847*bessely(0, (7940607665237127*r)/35184372088832))/2251799813685248, (747773142556017*besselk(0, (8879359236217549*r)/35184372088832))/1125899906842624 - (6083839523784019*besselj(0, (8879359236217549*r)/35184372088832))/18014398509481984 - (345*besseli(0, (8879359236217549*r)/35184372088832))/4503599627370496 + (6007480630199241*bessely(0, (8879359236217549*r)/35184372088832))/9007199254740992];

    diffarr = tarr-time;
    plot_time_index = find(diffarr == min(abs(diffarr)));
    if isempty(plot_time_index)
        plot_time_index = find(diffarr == -min(abs(diffarr)));
    end
    
    max_analytic = max(abs(feval(matlabFunction(annulus_equations(bessel_mode)),rarr)));
    max_numeric = max(abs(sol_matrix(plot_time_index,:)));
    
    if feval(matlabFunction(annulus_equations(bessel_mode)),rarr(2))*sol_matrix(plot_time_index,2) < 0
        scale_factor = -max_analytic/max_numeric;
    else
        scale_factor = -max_analytic/max_numeric;
    end
    
    figure(2);
    hold on;
    plot(rarr,scale_factor*sol_matrix(plot_time_index,:),'LineWidth',2);
    fplot(annulus_equations(bessel_mode),[Ri,Ro],'LineWidth',2);
    legend('Numeric Solution','Analytic Solution');
    xlabel('Radial Distance r (meters)');
    ylabel('Transverse Displacement u');
    hold off;
end

%% Algorithm Functions

%computes the next u(r,t_j+1) iteration
function nextu = computeNextu(r,u,h,high,k,beta,prevu,start_bc0_bool,start_bc0,start_bc1_bool,start_bc1,start_bc2_bool,start_bc2,start_bc3_bool,start_bc3,start_bc4_bool,start_bc4, ...
                                 end_bc0_bool,end_bc0,end_bc1_bool,end_bc1,end_bc2_bool,end_bc2,end_bc3_bool,end_bc3,end_bc4_bool,end_bc4)
                           
    laplapu = computelaplapu(r,u,h,high,start_bc1_bool,start_bc1,start_bc2_bool,start_bc2,start_bc3_bool,start_bc3,start_bc4_bool,start_bc4, ...
                             end_bc1_bool,end_bc1,end_bc2_bool,end_bc2,end_bc3_bool,end_bc3,end_bc4_bool,end_bc4);
                                  
    nextu = -(k^2)*laplapu/beta + 2*u - prevu;
    
    if start_bc0_bool
        nextu(1) = start_bc0;
    end
    if end_bc0_bool
        nextu(length(nextu)) = end_bc0;
    end
end

%computes the u(r,t) function for all times specified by the t array
function sol_matrix = computeChladniSolution(r,h,t,k,high,uinit,uinitinit,beta,start_bc0_bool,start_bc0,start_bc1_bool,start_bc1,start_bc2_bool,start_bc2,start_bc3_bool,start_bc3,start_bc4_bool,start_bc4, ...
                                             end_bc0_bool,end_bc0,end_bc1_bool,end_bc1,end_bc2_bool,end_bc2,end_bc3_bool,end_bc3,end_bc4_bool,end_bc4,save_to_csv,csv_filename,load_from_csv)
    % r - array of radial distance values to compute solution for
    % h - distance between distance values in r
    % t - array of times to compute solution for
    % k - distance between time values in t
    % high - 1 for sixth order accuracy finite differences, 0 for second
    % order accuracy finite differences
    % rinit - initial condition for the radial transverse displacement of
    % plate at t = 0
    % rinitinit - initial condition for the radial transverse displacement of
    % plate at t = -1
    % beta - material property
    % save_to_csv - true or false for saving the result into current MATLAB
    % directory
    % csv_filename - self explanatory
    % load_from_csv - determines whether the entire solution matrix is
    % cached in RAM (false) or not (true)
    % start/end_bcX_bool - indicates whether there is a boundary condition
    % at the inner or outer radius of order X. true or false
    % start/end_bcX - array of values (same length as t) for the boundary
    % condition at Ri or Ro of order X at each time.
    
    if ~load_from_csv
        sol_matrix = zeros(size(t,2),size(uinit,2));
        sol_matrix(1,:) = uinit;
    end
     
    if save_to_csv
        dlmwrite(csv_filename,uinit,'delimiter',',');
    end
    
    ucurr = uinit;
    uprev = uinitinit;
    
    for i = 1:length(t)-1
      
        sbc0 = start_bc0(i);
        sbc1 = start_bc1(i);
        sbc2 = start_bc2(i);
        sbc3 = start_bc3(i);
        sbc4 = start_bc4(i);
        ebc0 = end_bc0(i);
        ebc1 = end_bc1(i);
        ebc2 = end_bc2(i);
        ebc3 = end_bc3(i);
        ebc4 = end_bc4(i);
        
        %Print out the progress
        if mod(i,floor(length(t)/100)) == 0
            disp(sprintf('%3.2f%% done!',100*i/length(t)));
        end
        
        unext = computeNextu(r,ucurr,h,high,k,beta,uprev,start_bc0_bool,sbc0,start_bc1_bool,sbc1,start_bc2_bool,sbc2,start_bc3_bool,sbc3,start_bc4_bool,sbc4, ...
                                 end_bc0_bool,ebc0,end_bc1_bool,ebc1,end_bc2_bool,ebc2,end_bc3_bool,ebc3,end_bc4_bool,ebc4);            
                             
        if ~load_from_csv
            sol_matrix(i+1,:) = unext;
        end
        
        if save_to_csv
            dlmwrite(csv_filename,unext,'delimiter',',','-append');
        end
        
        uprev = ucurr;
        ucurr = unext;
    end 
end

%% Function for computing laplapu

%Computes the laplapu array from the u array as a function of radial radial
%distance. Can be either high precision (high = 1) or low precision (high = 0)
function laplapu = computelaplapu(r,u,h,high,start_bc1_bool,start_bc1,start_bc2_bool,start_bc2,start_bc3_bool,start_bc3,start_bc4_bool,start_bc4, ...
                       end_bc1_bool,end_bc1,end_bc2_bool,end_bc2,end_bc3_bool,end_bc3,end_bc4_bool,end_bc4)
    
    laplapu = zeros(1,length(u));

    for index = 1:length(u)
        curr_r = r(index);
        first = firstPartial(u,h,index,high,start_bc1_bool,start_bc1,end_bc1_bool,end_bc1);
        second = secondPartial(u,h,index,high,start_bc2_bool,start_bc2,end_bc2_bool,end_bc2);
        third = thirdPartial(u,h,index,high,start_bc3_bool,start_bc3,end_bc3_bool,end_bc3);
        fourth = fourthPartial(u,h,index,high,start_bc4_bool,start_bc4,end_bc4_bool,end_bc4);
        laplapu(index) = fourth + (2/curr_r)*third - (1/curr_r^2)*second + (1/curr_r^3)*first;
    end
end

%% Functions for computing partials
function res = firstPartial(u,h,index,high,bc_start_bool,bc_start,bc_end_bool,bc_end)
    if high == 1
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = ((3/4)*u(index+1)+(-3/20)*u(index+2)+(1/60)*u(index+3))/h;
            end
        elseif index == 2
            res = ((-3/4)*u(index-1)+(3/4)*u(index+1)+(-3/20)*u(index+2)+(1/60)*u(index+3))/h;
        elseif index == 3
            res = ((3/20)*u(index-2)+(-3/4)*u(index-1)+(3/4)*u(index+1)+(-3/20)*u(index+2)+(1/60)*u(index+3))/h;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = ((-1/60)*u(index-3)+(3/20)*u(index-2)+(-3/4)*u(index-1))/h;
            end
        elseif index == length(u)-1
            res =  ((-1/60)*u(index-3)+(3/20)*u(index-2)+(-3/4)*u(index-1)+(3/4)*u(index+1))/h;
        elseif index == length(u)-2
            res = ((-1/60)*u(index-3)+(3/20)*u(index-2)+(-3/4)*u(index-1)+(3/4)*u(index+1)+(-3/20)*u(index+2))/h;
        else
            res = ((-1/60)*u(index-3)+(3/20)*u(index-2)+(-3/4)*u(index-1)+(3/4)*u(index+1)+(-3/20)*u(index+2)+(1/60)*u(index+3))/h;
        end
    elseif high == 0
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = (1/2)*u(index+1)/h;
            end
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = (-1/2)*u(index-1)/h;
            end
        else
            res = ((-1/2)*u(index-1)+(1/2)*u(index+1))/h;
        end   
    end
end

function res = secondPartial(u,h,index,high,bc_start_bool,bc_start,bc_end_bool,bc_end)
    if high == 1
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = ((-49/18)*u(index)+(3/2)*u(index+1)+(-3/20)*u(index+2)+(1/90)*u(index+3))/h^2;
            end
        elseif index == 2
            res = ((3/2)*u(index-1)+(-49/18)*u(index)+(3/2)*u(index+1)+(-3/20)*u(index+2)+(1/90)*u(index+3))/h^2;
        elseif index == 3
            res = ((-3/20)*u(index-2)+(3/2)*u(index-1)+(-49/18)*u(index)+(3/2)*u(index+1)+(-3/20)*u(index+2)+(1/90)*u(index+3))/h^2;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = ((1/90)*u(index-3)+(-3/20)*u(index-2)+(3/2)*u(index-1)+(-49/18)*u(index))/h^2;
            end
        elseif index == length(u)-1
            res = ((1/90)*u(index-3)+(-3/20)*u(index-2)+(3/2)*u(index-1)+(-49/18)*u(index)+(3/2)*u(index+1))/h^2;
        elseif index == length(u)-2
            res = ((1/90)*u(index-3)+(-3/20)*u(index-2)+(3/2)*u(index-1)+(-49/18)*u(index)+(3/2)*u(index+1)+(-3/20)*u(index+2))/h^2;
        else
            res = ((1/90)*u(index-3)+(-3/20)*u(index-2)+(3/2)*u(index-1)+(-49/18)*u(index)+(3/2)*u(index+1)+(-3/20)*u(index+2)+(1/90)*u(index+3))/h^2;
        end
    elseif high == 0
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = (-2*u(index)+u(index+1))/h^2;
            end
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = (u(index-1)-2*u(index))/h^2;
            end
        else
            res = (u(index-1)-2*u(index)+u(index+1))/h^2;
        end   
    end
end

function res = thirdPartial(u,h,index,high,bc_start_bool,bc_start,bc_end_bool,bc_end)
    if high == 1
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = ((-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3)+(7/240)*u(index+4))/h^3;
            end
        elseif index == 2
            res = ((61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3)+(7/240)*u(index+4))/h^3;
        elseif index == 3
            res = ((-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3)+(7/240)*u(index+4))/h^3;
        elseif index == 4
            res = ((3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3)+(7/240)*u(index+4))/h^3;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = ((-7/240)*u(index-4)+(3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1))/h^3;
            end
        elseif index == length(u)-1
            res = ((-7/240)*u(index-4)+(3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1))/h^3;
        elseif index == length(u)-2
            res = ((-7/240)*u(index-4)+(3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2))/h^3;
        elseif index == length(u)-3
            res = ((-7/240)*u(index-4)+(3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3))/h^3;
        else
            res = ((-7/240)*u(index-4)+(3/10)*u(index-3)+(-169/120)*u(index-2)+(61/30)*u(index-1)+(-61/30)*u(index+1)+(-169/120)*u(index+2)+(-3/10)*u(index+3)+(7/240)*u(index+4))/h^3;
        end
    elseif high == 0
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = (-u(index+1)+(1/2)*u(index+2))/h^3;
            end
        elseif index == 2
            res = (u(index-1)-u(index+1)+(1/2)*u(index+2))/h^3;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = ((-1/2)*u(index-2)+u(index-1))/h^3;
            end
        elseif index == length(u)-1
            res = ((-1/2)*u(index-2)+u(index-1)-u(index+1))/h^3;
        else
            res = ((-1/2)*u(index-2)+u(index-1)-u(index+1)+(1/2)*u(index+2))/h^3;
        end   
    end
end

function res = fourthPartial(u,h,index,high,bc_start_bool,bc_start,bc_end_bool,bc_end)
    if high == 1
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = ((91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3)+(7/240)*u(index+4))/h^4;
            end
        elseif index == 2
            res = ((-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3)+(7/240)*u(index+4))/h^4;
        elseif index == 3
            res = ((169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3)+(7/240)*u(index+4))/h^4;
        elseif index == 4
            res = ((-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3)+(7/240)*u(index+4))/h^4;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = ((7/240)*u(index-4)+(-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index))/h^4;
            end
        elseif index == length(u)-1
            res = ((7/240)*u(index-4)+(-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1))/h^4;
        elseif index == length(u)-2
            res = ((7/240)*u(index-4)+(-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2))/h^4;
        elseif index == length(u)-3
            res = ((7/240)*u(index-4)+(-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3))/h^4;
        else
            res = ((7/240)*u(index-4)+(-2/5)*u(index-3)+(169/60)*u(index-2)+(-122/15)*u(index-1)+(91/8)*u(index)+(-122/15)*u(index+1)+(169/60)*u(index+2)+(-2/5)*u(index+3)+(7/240)*u(index+4))/h^4;
        end
    elseif high == 0
        if index == 1
            if bc_start_bool
                res = bc_start;
            else
                res = (6*u(index)-4*u(index+1)+u(index+2))/h^4;
            end
        elseif index == 2
            res = (-4*u(index-1)+6*u(index)-4*u(index+1)+u(index+2))/h^4;
        elseif index == length(u)
            if bc_end_bool
                res = bc_end;
            else
                res = (u(index-2)-4*u(index-1)+6*u(index))/h^4;
            end
        elseif index == length(u)-1
            res = (u(index-2)-4*u(index-1)+6*u(index)-4*u(index+1))/h^4;
        else
            res = (u(index-2)-4*u(index-1)+6*u(index)-4*u(index+1)+u(index+2))/h^4;
        end   
    end
end