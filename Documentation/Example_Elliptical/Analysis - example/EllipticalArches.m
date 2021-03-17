%% OpenArch - Script to plot an admissible thrust line of an elliptical arch

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
t_over_R = 0.20602; %calculate this value using other script
Thickness = Rise*t_over_R;
Num_Blocks = 20;
RuptureAngles = [NaN,NaN,NaN]; %Take from analytical results
ReactionLocations = [-Rise/(2*c),Rise/(2*c)]; %x coordinates of left and right reactions


% Input Variables
density = 1;     
L = 0;          %amount of additional mass (L = 0: no additional mass, L = 1: 1x mass of the arch)
alpha_m = 0.1;  %inertial contribution of the additional mass
%%
%--------------------------------------------------------------------------
%---------------------------- G E O M E T R Y -----------------------------
%--------------------------------------------------------------------------

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

[xin,yin] = Toolkit.pol2car(intrados,thetas_int);
[xex,yex] = Toolkit.pol2car(extrados,thetas_ext);

mid  = (intrados+extrados)/2;
[xmid,ymid] = Toolkit.pol2car(mid,thetas_int);


% figure
% hold on
% axis equal
% axis off
% 
% plot(xmid,ymid,'k--','linewidth',1)
% plot(xin,yin,'k','linewidth',2)
% plot(xex,yex,'k','linewidth',2)
% plot([xin(end) xex(end)],[yin(end) yex(end)],'k','linewidth',2)
% plot([xin(1) xex(1)],[yin(1) yex(1)],'k','linewidth',2)
% hold off

%%
%--------------------------------------------------------------------------
%---------------------- L A T E R A L   L O A D I N G ---------------------
%--------------------------------------------------------------------------

figure
%Estimated reactions from virtual work
[Rv_L,Rh_L,Rv_R,Rh_R,M,thetas] = Toolkit.PinPinReactions(ARCH);

[x_ext,y_ext] = Toolkit.pol2car(ARCH(1).r(3),ARCH(1).theta(3));
[x_int,y_int] = Toolkit.pol2car(ARCH(1).r(1),ARCH(1).theta(1));
t_support = x_ext-x_int;

%right hand side starting point - specify a range or a single value
% StartingPoints = linspace(x_ext,x_ext-1*t_support,15);
StartingPoints = x_ext;

for k = 1 :length(StartingPoints)
    StartingPoint = [StartingPoints(k),0];
    
    % form load line
    global xp; % x coordinate of pole
    global yp; % y coordinate of pole
    global Force_Polygon; %plot of the force polygon
    global Thrust_Line; %plot of arch and thrust line
    
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
    Force_Polygon = subplot(1,4,4);
    Toolkit.ForcePolygon(x_ll,y_ll)
    
    
    Thrust_Line = subplot(1,4,1:3);
    Toolkit.Build(ARCH)
    hold on
    [x,y] = Toolkit.ThrustLine(ARCH,x_ll,y_ll);
    
    %%
    %--------------------------------------------------------------------------
    %----------------- M I N I M U M   T H R U S T   L I N E ------------------
    %--------------------------------------------------------------------------
    [x_ext,~] = Toolkit.pol2car(ARCH(end).r(4),ARCH(end).theta(4));
    [x_int,~] = Toolkit.pol2car(ARCH(end).r(2),ARCH(end).theta(2));
    
    %left hand side end point - specify a range or a single value
    xfinds = linspace(x_ext,x_int,200);
%     xfinds = x_int;
%          xfinds = x_ext;
    for i = 1:length(xfinds)
        
        
        tol = 0.0000000001;
        [x,y,success] = PoleFinding.find_minimum_thrust(ARCH,x,y,x_ll,y_ll,xfinds(i),tol);
        if success == "yes"
            break
        end
    end
    
    if success == "yes"
        delete(Force_Polygon);
        delete(Thrust_Line);
        Force_Polygon = subplot(1,4,4);
        
        Toolkit.ForcePolygon(x_ll,y_ll)
        
        Thrust_Line = subplot(1,4,1:3);
        Toolkit.Build(ARCH)
        hold on
        plot(x,y,'b-','linewidth',2)
        break
    end
    delete(Force_Polygon);
    delete(Thrust_Line);
end

if success == "no"
    error('INCREASE THICKNESS')
end


%%
%--------------------------------------------------------------------------
%---------------------- P L O T   R E S U L T S ---------------------------
%--------------------------------------------------------------------------

figure
hold on
axis equal
axis off

% plot arch and thrust line:
% x = linspace(Rise/(2*c),-Rise/(2*c),Num_Blocks+1); 
% y = Catenary(Rise,x);
plot(xmid,ymid,'k--','linewidth',1)
plot(xin,yin,'k','linewidth',2)
plot(xex,yex,'k','linewidth',2)
plot([xin(end) xex(end)],[yin(end) yex(end)],'k','linewidth',2)
plot([xin(1) xex(1)],[yin(1) yex(1)],'k','linewidth',2)
plot(x,y,'b','linewidth',2)

%plot hinge locations:
hinges = Toolkit.HingeLocations(ARCH,x,y,5,"extrados");

for i = 1:length(hinges)
    xh = x(hinges(i));
    yh = y(hinges(i));
    plot(xh,yh,'k.','markersize',14);
    plot([0,xh],[0,yh],'k--','linewidth',1)
end

hold off

%% Save results to folder
Result = struct;
Result.geom = ARCH;
Result.x = x;
Result.y = y;

%Give file a sensible name e.g. Circ03 - circular arch subject to 0.3*9.81
%lateral acceleration
save(['Circ03.mat'],'Result');