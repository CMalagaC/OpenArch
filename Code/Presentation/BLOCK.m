%% OpenArch - BLOCK Class definition

%  Contact:
%  T. McLean, thomas.o.mclean@outlook.com
%  C. Malaga-Chuquitaype, c.malaga@imperial.com

%% Notes

%  The coordinates 'r' and 'theta' should be entered such that the first 2
%  sets of coordinates correspond to the intrados of the arch, and the last
%  2 sets correspond to the extrados in order of increasing 'theta'.

%%
classdef BLOCK
    properties
        
        %input properties
        r;               %[m] radial coordinates of 4 corners of trapezoid
        theta;           %[theta] angular coordinates of 4 corners of trapezoid
        density;         %[kg/m3] material density
        
        %calculated properties
        A;               %[m2] area of block
        M;               %[kg] mass of block
        W;               %[N] weight of block
        
        r_CG;            %[m] radial coordinate of centre of gravity
        theta_CG;        %[rad] angular coordinate of centre of gravity
        
        Fv;              % External vertical force applied at centroid (downwards positive)
        Fh               % External horizontal force applied at centroid (left to right positive)
        
    end
    
    methods
        
        function obj = BLOCK(r,theta,density)
            
            global GravitationalAcceleration
            global LateralAcceleration
            
            %% Create object
            obj.r = r;
            obj.theta = theta;
            obj.density = density;
            
            %% Calculate area of the block and centre of gravity
            %Convert coordinates to cartesian for simplicity
            [x,y] = Toolkit.pol2car(r,theta);
            %reorder x and y
            x = [x(1) x(2) x(4) x(3)];
            y = [y(1) y(2) y(4) y(3)];
            
            A = 0;
            x_cg = 0;
            y_cg = 0;
            for i = 1:length(x)-1
                A = A + (x(i)*y(i+1)-x(i+1)*y(i));
                x_cg = x_cg + (x(i)+x(i+1))*(x(i)*y(i+1)-x(i+1)*y(i));
                y_cg = y_cg + (y(i)+y(i+1))*(x(i)*y(i+1)-x(i+1)*y(i));
            end

            obj.A = abs(A/2);
            
            x_cg = -x_cg/(6*obj.A);
            y_cg = -y_cg/(6*obj.A);
            
            [obj.r_CG,obj.theta_CG]=Toolkit.car2pol(x_cg,y_cg);
            
            %% Calculate Mass of the block
            obj.M = obj.A*density;
            
            %% Calculate Weight of the block
            obj.W = obj.M*GravitationalAcceleration;
            
            %% Calculate Applied horizontal force (equivalent seismic loading)
            obj.Fh = LateralAcceleration*obj.M; 
            obj.Fv = 0; %set to zero initially;
            
        end
        
        
    end
end


