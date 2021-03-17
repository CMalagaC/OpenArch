%% OpenArch - Script to find the minimum thickness of an elliptical arch

%  Contact:
%  T. McLean, thomas.o.mclean@outlook.com
%  C. Malaga-Chuquitaype, c.malaga@imperial.com

clear
clc
close all


%%
%--------------------------------------------------------------------------
%--------------------------- V A R I A B L E S ----------------------------
%--------------------------------------------------------------------------

%Global Variables
global GravitationalAcceleration;
global LateralAcceleration;
global Rise;
global c; %rise/span
global Thickness;
global Num_Blocks;
global RuptureAngles;
global ReactionLocations
global StartingPoint

GravitationalAcceleration = 1*9.81;
LateralAcceleration = 0.3*9.81;
Rise = 10;
c = 1/2;
t_over_R = 0.5; %input an overestimate here
Thickness = Rise*t_over_R;
Num_Blocks = 20;
RuptureAngles = [NaN,NaN,NaN]; %Take from analytical results
ReactionLocations = [-Rise/(2*c),Rise/(2*c)]; %x coordinates of left and right reactions


% Input Variables
density = 1;     
L = 0;          %amount of additional mass (L = 0: no additional mass, L = 1: 1x mass of the arch)
alpha_m = 0;  %inertial contribution of the additional mass

dts = [0.1, 0.01, 0.001, 0.0001, 0.00001]; %increments to reduce thickness by
%%
%--------------------------------------------------------------------------
%---------------------------- G E O M E T R Y -----------------------------
%--------------------------------------------------------------------------
tstart = tic;
for m = 1:length(dts)
    dt = dts(m);
    success = "yes";
    while success == "yes"
        t_over_R = t_over_R - dt;
        if t_over_R < 0.00001
            break
        end
        Thickness = Rise*t_over_R;
        
        b = Rise;
        a = b/(2*c);
        
        %Generate elliptical geometry:
        %intrados
        a = a-Thickness/2;
        b = b-Thickness/2;
        thetas_int = linspace(0,pi,Num_Blocks+1);
        intrados = (a*b)./sqrt((b*cos(thetas_int)).^2+(a*sin(thetas_int)).^2);
        
        
        b = Rise;
        a = b/(2*c);
        %extrados
        a = a+Thickness/2;
        b = b+Thickness/2;
        thetas_ext = linspace(0,pi,Num_Blocks+1);
        extrados = (a*b)./sqrt((b*cos(thetas_ext)).^2+(a*sin(thetas_ext)).^2);
        
        % Assign coordinates to BLOCKs
        ARCH=BLOCK.empty;
        
        index = [0,1];
        for i = 1:Num_Blocks
            index = index+1;
            
            ARCH(end+1) = BLOCK([intrados(index),extrados(index)],[thetas_int(index),thetas_ext(index)],density);
        end
        
        % add soil mass
        ARCH = AddSoil(ARCH,L,alpha_m);
        
        
        %%
        %--------------------------------------------------------------------------
        %---------------------- L A T E R A L   L O A D I N G ---------------------
        %--------------------------------------------------------------------------
        
        %Estimated reactions from vertical work
        [Rv_L,Rh_L,Rv_R,Rh_R,M,thetas] = Toolkit.PinPinReactions(ARCH);
       

        [x_ext,y_ext] = Toolkit.pol2car(ARCH(1).r(3),ARCH(1).theta(3));
        [x_int,y_int] = Toolkit.pol2car(ARCH(1).r(1),ARCH(1).theta(1));
        t_support = x_ext-x_int;
        
        %right hand side starting point - specify a range or a single value
%             StartingPoints = linspace(x_ext,x_ext-1*t_support,200);
      StartingPoints = x_ext;
        
        for k = 1 :length(StartingPoints)
            StartingPoint = [StartingPoints(k),0];
            
            % form load line
            global xp; % x coordinate of pole
            global yp; % y coordinate of pole
            
            %Create load line
            %use reactions from virtual work to calculate where to start the load line
            x_ll = zeros(1,Num_Blocks+1);
            y_ll = zeros(1,Num_Blocks+1);
            x_ll(1) = Rh_L;
            y_ll(1) = Rv_L;
            
            for i = 1:Num_Blocks
                x_ll(i+1) = x_ll(i) + ARCH(i).Fh;
                y_ll(i+1) = y_ll(i) - (ARCH(i).W+ARCH(i).Fv);
            end
            %reorder load line from right reaction to left reaction
            x_ll = x_ll(end:-1:1);
            y_ll = y_ll(end:-1:1);
            
            
            %Set pole at (0,0) in load space (arbitrary);
            xp = -5000; yp = -5000; %for high levels of loading place at (-5000,-5000) or more to ensure convergence
            
            [x,y] = Toolkit.ThrustLine(ARCH,x_ll,y_ll);
            
            %%
            %--------------------------------------------------------------------------
            %----------------- M I N I M U M   T H R U S T   L I N E ------------------
            %--------------------------------------------------------------------------
            [x_ext,~] = Toolkit.pol2car(ARCH(end).r(4),ARCH(end).theta(4));
            [x_int,~] = Toolkit.pol2car(ARCH(end).r(2),ARCH(end).theta(2));
            t = abs(x_ext-x_int); %use this t to specify required range of values within support
            
          %left hand side end point - specify a range or a single value if
          %known
%         xfinds = linspace(x_int,x_ext,200);
          xfinds = linspace(x_ext,x_int,200);
%          xfinds = x_ext;
            
            for i = 1:length(xfinds)
                tol = 0.0000000001;
                [x,y,success] = PoleFinding.find_minimum_thrust(ARCH,x,y,x_ll,y_ll,xfinds(i),tol);
                if success == "yes"
                    break
                end
            end
            
            if success == "yes"
                break
            end
        end
        
    end
    t_over_R = t_over_R + dt;
    
end
timing = toc(tstart);
fprintf('Min Thickness = %g\n',t_over_R)
fprintf('Time taken = %.3f\n',timing)