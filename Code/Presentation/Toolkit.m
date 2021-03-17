%% OpenArch - Compilation of useful functions

%  Contact:
%  T. McLean, thomas.o.mclean@outlook.com
%  C. Malaga-Chuquitaype, c.malaga@imperial.com

classdef Toolkit
    
    methods (Static)
        %% GENERAL
        %Polar coordinates to cartesian
        function [x,y] = pol2car(r,theta)
            x = r.*cos(theta);
            y = r.*sin(theta);
        end
        
        %Cartesian coordinates to polar
        function [r,theta] = car2pol(x,y)
            r = sqrt(x.^2 + y.^2);
            theta = atan2(y,x);
        end
        
        %Creates an arch of constant thickness given a midline (x,y)
        function [intrados,extrados,thetas] = ConstantThickness(x,y,Thickness)
            
            intrados = zeros(size(x));
            extrados = zeros(size(x));
            thetas = zeros(size(x));
            
            for i = 2:length(x)-1
                [r,theta] = Toolkit.car2pol(x(i),y(i));
                m1 = (y(i)-y(i+1))/(x(i)-x(i+1));
                m2 = (y(i)-y(i-1))/(x(i)-x(i-1));
                m = (m1+m2)/2;
                
                theta_mid = atan2(m,1);
                phi = pi/2 - theta + theta_mid;
                delta = Thickness/cos(phi);
                
                intrados(i) = r - delta/2;
                extrados(i) = r + delta/2;
                thetas(i) = theta;
            end
            
            mr = (y(2)-y(1))/(x(2)-x(1));
            ml = (y(end)-y(end-1))/(x(end)-x(end-1));
            
            %right hand side
            [x,y] = Toolkit.pol2car(extrados(2),thetas(2));
            c_ex = y-mr*x;
            y_ex = 0; x_ex = -c_ex/mr;
            [extrados(1),thetas(1)] = Toolkit.car2pol(x_ex,y_ex);
            
            [x,y] = Toolkit.pol2car(intrados(2),thetas(2));
            c_in = y-mr*x;
            y_in = 0; x_in = -c_in/mr;
            [intrados(1),~] = Toolkit.car2pol(x_in,y_in);
            
            %left hand side
            [x,y] = Toolkit.pol2car(extrados(end-1),thetas(end-1));
            c_ex = y-ml*x;
            y_ex = 0; x_ex = -c_ex/ml;
            [extrados(end),thetas(end)] = Toolkit.car2pol(x_ex,y_ex);
            
            [x,y] = Toolkit.pol2car(intrados(end-1),thetas(end-1));
            c_in = y-ml*x;
            y_in = 0; x_in = -c_in/ml;
            [intrados(end),~] = Toolkit.car2pol(x_in,y_in);
            
        end
        
        %Plot blocks and their centres of gravity
        function Build(ARCH)
            global Thickness;
            global Rise;
            global LateralAcceleration
            
            hold on
            grid on
            axis equal tight
            %             xlim([-1.5 1.5])
            %             ylim([0 1.5])
            for i = 1:length(ARCH)
                thetas_int = ARCH(i).theta(1:2);
                thetas_ext = ARCH(i).theta(3:4);
                r_int = ARCH(i).r(1:2);
                r_ext = ARCH(i).r(3:4);
                r_CG = ARCH(i).r_CG;
                theta_CG = ARCH(i).theta_CG;
                
                [x_int,y_int] = Toolkit.pol2car(r_int,thetas_int);
                [x_ext,y_ext] = Toolkit.pol2car(r_ext,thetas_ext);
                
                %Fill in block
                fill([x_int(1),x_int(2),x_ext(2),x_ext(1)],[y_int(1),y_int(2),y_ext(2),y_ext(1)],'k')
                alpha(0.1)
                
                %                 plot intrados
                %                                 plot(x_int,y_int,'-k','linewidth',2);
                
                %                 plot extrados
                %                                 plot(x_ext,y_ext,'-k','linewidth',2);
                
                %plot centre of mass
                [x_CG,y_CG] = Toolkit.pol2car(r_CG,theta_CG);
                %                 plot(x_CG,y_CG,'k.','markersize',12)
                
                %
                %drawnow
                %pause(0.1)
            end
            xlabel('x [m]')
            ylabel('y [m]')
            %             title([num2str(LateralAcceleration/9.81),'g,  t/R = ',num2str(Thickness/Rise)])
            %             h(1) = plot(NaN,NaN,'k.','markersize',6);
            %             legend(h,'Centre of gravity')
            %             legend boxoff
            
            
            
        end
        
        %Function to plot force polygon with current pole position (xp,yp)
        function ForcePolygon(x_ll,y_ll)
            global xp
            global yp
            hold on
            axis equal
            
            grid on
            %             grid minor
            plot(xp,yp,'.r','markersize',40)
            plot(x_ll,y_ll,'-*k','linewidth',2)
            for i = 1:length(x_ll)
                plot([xp,x_ll(i)],[yp,y_ll(i)],'b','linewidth',2)
            end
            set(gca,'fontsize',15)
            hold off
        end
        
        %Total vertical load on the arch (weight + external)
        % NOTE: for this function downards is defined as positive
        function W = TotalVertical(ARCH)
            W = 0;
            for i = 1:length(ARCH)
                W = W + ARCH(i).W + ARCH(i).Fv;
            end
        end
        
        %Total horizontal load on the arch (external)
        function H = TotalHorizontal(ARCH)
            H = 0;
            for i = 1:length(ARCH)
                H = H + ARCH(i).Fh;
            end
        end
        
        %Total applied moment
        function M_total=AppliedMoment(ARCH)
            global Num_Blocks
            global ReactionLocations;
            
            % Take moment about right side
            M_total = 0;
            for i = 1:Num_Blocks
                W = ARCH(i).W;
                Fv = ARCH(i).Fv;
                Fh = ARCH(i).Fh;
                r = ARCH(i).r_CG;
                theta = ARCH(i).theta_CG;
                
                M_total = M_total + (W+Fv)*(ReactionLocations(2)-r*cos(theta)) - Fh*r*sin(theta);
            end
        end
        
        %% REACTIONS
        %Calculate reactions of pinned structure due to horizontal and
        %vertical loading (no applied moment)
        
        function [Rv_L,Rh_L,Rv_R,Rh_R,M,thetas] = PinPinReactions(ARCH)
            
            global Num_Blocks;
            global ReactionLocations;
            
            
            %% Calculate right side vertical reaction through moment equilibrium
            
            % Take moment about right side
            moment = 0;
            for i = 1:Num_Blocks
                W = ARCH(i).W;
                Fv = ARCH(i).Fv;
                Fh = ARCH(i).Fh;
                r = ARCH(i).r_CG;
                theta = ARCH(i).theta_CG;
                
                moment = moment + (W+Fv)*(ReactionLocations(2)-r*cos(theta)) - Fh*r*sin(theta);
            end
            
            % Vertical equilibrium
            Rv_L = moment/(ReactionLocations(2)-ReactionLocations(1));
            Rv_R = Toolkit.TotalVertical(ARCH) - Rv_L;
            
            %% Principle of virtual work for right hand reaction
            %             W_total = Toolkit.TotalVertical(ARCH);
            %             Rv = W_total/2;
            
            M = zeros(1,length(ARCH)+2);
            m = zeros(1,length(ARCH)+2);
            thetas = zeros(1,length(ARCH)+2); thetas(end) = pi;
            
            for i = 1:Num_Blocks
                %centre of gravity of block
                r_i = ARCH(i).r_CG;
                theta_i = ARCH(i).theta_CG;
                
                M(i+1) = Rv_R*(ReactionLocations(2) -r_i*cos(theta_i));
                
                % Moments due to loading on previous blocks
                if i ~=1
                    for j = 1:i-1
                        %centre of gravity of previous block
                        r_j = ARCH(j).r_CG;
                        theta_j = ARCH(j).theta_CG;
                        %weight + external load of previous block
                        Wj = ARCH(j).W + ARCH(j).Fv;
                        Fh = ARCH(j).Fh;
                        
                        M(i+1) = M(i+1) -Wj*(r_j*cos(theta_j)-r_i*cos(theta_i)) + Fh*(r_i*sin(theta_i)-r_j*sin(theta_j));
                        
                    end
                end
                
                m(i+1) = r_i*sin(theta_i); %moment due to unit horizontal load at right support (right to left)
                thetas(i+1) = theta_i;
            end
            
            %Calculate integral I = int(M(theta)*m(theta))dtheta
            I = trapz(thetas,M.*m);
            
            %Calculate integral II = int(m(theta)*m(theta))dtheta
            II = trapz(thetas,m.*m);
            
            X = -I/II; %compatability condition
            Rh_R = 1*X;
            Rh_L = -Toolkit.TotalHorizontal(ARCH) - Rh_R;
            
            M = M + X*m;
            
        end
        
        %% THRUST LINES
        %Calculate thrust line coordinates using rotational equilibrium
        function [x,y] = ThrustLine(ARCH,x_ll,y_ll)
            global xp
            global yp
            global StartingPoint
            
            %Calculate thrust forces from the force polygon
            Th = zeros(1,length(x_ll));
            Tv = zeros(1,length(x_ll));
            
            for i = 1:length(x_ll)
                Th(i) = x_ll(i)-xp;
                Tv(i) = yp-y_ll(i);
            end
            
            %Initialise thrust line coordinates
            x = zeros(1,length(ARCH)+1);
            y = zeros(1,length(ARCH)+1);
            x(1) = StartingPoint(1);
            y(1) = StartingPoint(2);
            
            %Calculate thrust line coordinates
            for i = 1:length(ARCH)
                block = ARCH(i);
                
                x1 = x(i); y1 = y(i);
                [xcg,ycg] = Toolkit.pol2car(block.r_CG,block.theta_CG);
                Th1 = Th(i); Tv1 = Tv(i);
                Th2 = Th(i+1); Tv2 = Tv(i+1);
                
                %equation of voussoir edge
                [xin,yin] = Toolkit.pol2car(block.r(2),block.theta(2));
                [xex,yex] = Toolkit.pol2car(block.r(4),block.theta(4));
                m = (yex-yin)/(xex-xin);
                c = yin-m*xin;
                
                %formulate matrix equations
                a = Tv2*xcg + Th2*ycg + Tv1*(x1-xcg) - Th1*(ycg-y1);
                
                A = [a; -c];
                T = [Tv2 Th2; m -1];
                
                X = T\A;
                
                x(i+1) = X(1);
                y(i+1) = X(2);
                
            end
            
        end
        
        %Check if the thrust line is with the arch geometry
        function in_or_out = CheckThrust(ARCH,x,y)
            % in_or_out returns:
            %                   - "in" if thrust line inside extrados and
            %                   intrados
            %                   - "out" if thrust line outside extrados and
            %                   intrados
            
            decisions = zeros(size(x));
            r = Toolkit.car2pol(x,y);
            for i = 2:length(x)-1
                r_int = ARCH(i-1).r(2);
                r_ext = ARCH(i-1).r(4);
                
                if r(i) - r_int< -0.00000000000001 || r(i)-r_ext > 0.00000000000001
                    decisions(i) = 1;
                end
            end
            
            if sum(decisions) == 0  %whole line inside
                in_or_out = "in";
            else
                in_or_out = "out";
            end
            
        end
        
        function target_height = Target(TL_coords)
            %find target height for SW line
            for i = 1:length(TL_coords.x)-1
                if abs(TL_coords.x(i)) < 0.0001
                    target_height = TL_coords.y(i);
                    break
                elseif TL_coords.x(i)*TL_coords.x(i+1)<1
                    x1 = TL_coords.x(i); y1 = TL_coords.y(i);
                    x2 = TL_coords.x(i+1); y2 = TL_coords.y(i+1);
                    
                    m = (y2-y1)/(x2-x1);
                    c = y1 - m*x1;
                    
                    target_height = c;
                    
                    break
                end
            end
        end
        
        %determines whether the thrust line lies above or below the arch
        function location = AboveOrBelow(Geom,x,y)
            
            
            decisions = zeros(size(x)); %1 = above, 0 = inside, -1 = below
            errors = zeros(size(x)); %absolute vertical distance
            
            for i = 2:length(x)-1
                block = Geom(i-1);
                [~,y1] = Toolkit.pol2car(block.r(2),block.theta(2));
                [~,y2] = Toolkit.pol2car(block.r(4),block.theta(4));
                y_max = max([y1,y2]);
                y_min = min([y1,y2]);
                if y(i) >y_max
                    decisions(i) = 1;
                    errors(i) = y(i)-y_max;
                elseif y(i) < y_min
                    decisions(i) = -1;
                    errors(i) = y_min-y(i);
                end
            end
            
            %find position of biggest difference from geometry
            index = find(errors == max(errors));
            if decisions(index) == 1
                location = "above";
            elseif decisions(index) == -1
                location = "below";
            end
            
            
        end
        
        %Minimum distance of thrust line to arch extremes (for admissible
        %thrust lines)
        function Min_diff = min_of_diffs(ARCH,x,y)
            diffs = zeros(size(x));
            for i = 2:length(x)-1
                %find coordinates of current voussoir edge for check
                [x_int,y_int] = Toolkit.pol2car(ARCH(i-1).r(2),ARCH(i-1).theta(2));
                [x_ext,y_ext] = Toolkit.pol2car(ARCH(i-1).r(4),ARCH(i-1).theta(4));
                
                diff_int = sqrt((x(i)-x_int)^2+(y(i)-y_int)^2);
                diff_ext = sqrt((x(i)-x_ext)^2+(y(i)-y_ext)^2);
                diffs(i) = min(diff_int,diff_ext);
                
            end
            
            Min_diff = min(diffs(2:end-1));
        end
        
        %Maximum distance of thrust line to arch extremes (for non admissible
        %thrust lines)
        function Max_diff = max_of_diffs(ARCH,x,y)
            diffs = zeros(size(x));
            for i = 2:length(x)-1
                %find coordinates of current voussoir edge for check
                [x_int,y_int] = Toolkit.pol2car(ARCH(i-1).r(2),ARCH(i-1).theta(2));
                [x_ext,y_ext] = Toolkit.pol2car(ARCH(i-1).r(4),ARCH(i-1).theta(4));
                
                
                x_min = min(x_int,x_ext);
                x_max = max(x_int,x_ext);
                
                y_min = min(y_int,y_ext);
                y_max = max(y_int,y_ext);
                
                if x(i)>x_max || x(i)<x_min
                    if y(i)>y_max || y(i)<y_min
                        diff_int = sqrt((x(i)-x_int)^2+(y(i)-y_int)^2);
                        diff_ext = sqrt((x(i)-x_ext)^2+(y(i)-y_ext)^2);
                        diffs(i) = min(diff_int,diff_ext); %outside intrados and extrados
                    end
                end
                
                
            end
            
            Max_diff = max(diffs);
        end
        
        %Function to find the hinge locations given an ARCH geometry and a
        %limit thrust line
        function [hinge_locations] = HingeLocations(ARCH,x,y,n,start)
            
            %n = number of hinges to look for
            % start = "intrados" first hinge at intrados
            % start = "extrados" first hinge at extrados
            %% get intrados and extrados
            
            intrados = zeros(1,length(ARCH)+1);
            extrados = zeros(1,length(ARCH)+1);
            thetas_int = zeros(1,length(ARCH)+1);
            thetas_ext = zeros(1,length(ARCH)+1);
            
            for i =1:length(ARCH)
                intrados(i) = ARCH(i).r(1);
                extrados(i) = ARCH(i).r(3);
                thetas_int(i) = ARCH(i).theta(1);
                thetas_ext(i) = ARCH(i).theta(3);
            end
            intrados(end) = ARCH(end).r(2);
            extrados(end) = ARCH(end).r(4);
            thetas_int(end) = ARCH(end).theta(2);
            thetas_ext(end) = ARCH(end).theta(4);
            
            [xin,yin] = Toolkit.pol2car(intrados,thetas_int);
            
            %intrados modified
            i = 1;
            while yin(i) == yin(i+1)
                yin(i) = NaN;
                xin(i) = NaN;
                yin(end-(i-1)) = NaN;
                xin(end-(i-1)) = NaN;
                i = i+1;
            end
            
            [xex,yex] = Toolkit.pol2car(extrados,thetas_ext);
            
            %% Calculate distances from thrust line to intrados and extrados
            diffs_ex = zeros(size(x));
            diffs_in = zeros(size(x));
            
            for i = 1:length(x)
                diffs_in(i) = sqrt((xin(i)-x(i))^2+(yin(i)-y(i))^2);
                diffs_ex(i) = sqrt((xex(i)-x(i))^2+(yex(i)-y(i))^2);
                
            end
            
            
            
            hinge_locations = zeros(1,n);
            
            %Determine where to start looking (goes from right to left)
            if start == "extrados"
                looking = "extrados";
            elseif start == "intrados"
                looking = "intrados";
            end
            
            Limit_index = 1;
            for i = 1:n % for each hinge
                
                if looking == "extrados"
                    for j = Limit_index:length(x)-1
                        if diffs_ex(j+1)>diffs_ex(j)
                            if diffs_ex(j)~= NaN
                                if diffs_ex(j)<0.1
                                    Limit_index = j+1;
                                    break
                                end
                            end
                            
                        end
                    end
                    hinge_locations(i) = Limit_index-1;
                    looking = "intrados";
                    
                    
                elseif looking == "intrados"
                    for j = Limit_index:length(x)-1
                        if diffs_in(j+1)>diffs_in(j)
                            if diffs_in(j)~= NaN
                                if diffs_in(j)<0.1
                                    Limit_index = j+1;
                                    break
                                end
                            end
                            
                        end
                    end
                    hinge_locations(i) = Limit_index-1;
                    looking = "extrados";
                end
                
                %check extrados and intrados hinge at left hand side
                
                if diffs_ex(end)< 0.0001
                    hinge_locations(end) = length(x);
                end
                
                if diffs_in(end)<0.0001
                    hinge_locations(end) = length(x);
                end
                
            end
            
            
            
        end
        
        
        
    end
    
end



