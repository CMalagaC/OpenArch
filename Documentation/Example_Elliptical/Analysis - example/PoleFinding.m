%% Individual Research Project
% Thomas McLean CID: 01198899
% Functions for finding position of the pole of force polygon

classdef PoleFinding
    
    methods (Static)
        
        function [x,y] = Move_xy(dx,dy,x_ll,y_ll,ARCH)
            global xp
            global yp
         
            
            xp = xp + dx;
            yp = yp + dy;
            
            
            
            [x,y] = Toolkit.ThrustLine(ARCH,x_ll,y_ll);

            
        end
        
        function [x,y] = find_end_thrust(ARCH,x,y,x_ll,y_ll,tol,x_find)
            %updates global xp and yp (position of pole in 'force space')
            %such that the thrust line passes through the coordinate x_left
            global xp;
            global yp;
            
            
            %take initial increment as 10% of left vertical reaction
            dy = 0.10*abs(y_ll(1));
            
            
            if x(end) < x_find %thrust line to the left of the intrados
                
                %move pole until the line starts to the right of the
                %intrados
                while x(end) < x_find
                    [x,y] = PoleFinding.Move_xy(0,-dy,x_ll,y_ll,ARCH);
                end
                yp_right = yp; %pole that gives thrust line right of intrados
                yp_left = yp+dy;
                
                error = abs(x_find-x(end));
                while error>tol
                    %use bisection method to iterate to accurate solution
                    
                    yp_new = (yp_left+yp_right)/2;
                    dy = yp_new - yp;
                    [x,y] = PoleFinding.Move_xy(0,dy,x_ll,y_ll,ARCH);
                    
                    if abs(x(end)-x_find)<0.00001
                        %                         [x,y] = PoleFinding.Move_xy(0,0,x_ll,y_ll,ARCH);
                        yp=yp;
                        xp=xp;
                        return
                    end
                    
                    %update either yp_left or yp_right
                    if x(end)<x_find
                        yp_left = yp;
                    elseif x(end) >x_find
                        yp_right = yp;
                    end
                    error = abs(x_find-x(1));
                end
                
                %make sure thrust line starts inside intrados:
                yp_new = yp_left;
                dy = yp_new - yp;
                [x,y] = PoleFinding.Move_xy(0,dy,x_ll,y_ll,ARCH);
                
                
            elseif x(end) >x_find %thrust line to the right of the intrados
                
                %move pole until the line starts to the left of the
                %intrados
                while x(end) > x_find
                    [x,y] = PoleFinding.Move_xy(0,dy,x_ll,y_ll,ARCH);
                end
                yp_left = yp; %pole that gives thrust line right of intrados
                yp_right = yp-dy;
                
                
                error = abs(x_find-x(end));
                while error>tol
                    %use bisection method to iterate to accurate solution
                    
                    yp_new = (yp_left+yp_right)/2;
                    dy = yp_new - yp;
                    [x,y] = PoleFinding.Move_xy(0,dy,x_ll,y_ll,ARCH);
                    
                    if abs(x(end)-x_find)<0.00001
                        [x,y] = PoleFinding.Move_xy(0,0,x_ll,y_ll,ARCH);
                        return
                    end
                    
                    %update either yp_left or yp_right
                    if x(end)<x_find
                        yp_left = yp;
                    elseif x(end) >x_find
                        yp_right = yp;
                    end
                    error = abs(x_find-x(end));
                end
                
                %make sure thrust line starts inside intrados:
                yp_new = yp_left;
                dy = yp_new - yp;
                [x,y] = PoleFinding.Move_xy(0,dy,x_ll,y_ll,ARCH);
            end
            
            
        end
        
        function [x,y,success,error] = find_minimum_thrust(ARCH,x,y,x_ll,y_ll,x_find,tol)
            %Updates position of pole to find minimum thrust line
            global xp;
            global yp;
            
            
            %move thrust line so that it lies at left intrados
            [x,y] = PoleFinding.find_end_thrust(ARCH,x,y,x_ll,y_ll,tol,x_find);
            
            %move pole until the line is within the intrados and extrados
            dx = max(0.1*abs(x_ll(1)),0.05*abs(y_ll(1)));
            
            
            check = Toolkit.CheckThrust(ARCH,x,y);
            if check == "out"
                error = Toolkit.max_of_diffs(ARCH,x,y);
            elseif check == "in"
                error = Toolkit.min_of_diffs(ARCH,x,y);
            end
            
            dxs = 0;
            
            if check == "out"
                
                %move pole to the left or right depending on position
                %of the thrust line (increase of decrease horizontal
                %thrust)
                start_location = Toolkit.AboveOrBelow(ARCH,x,y);
                while check == "out"
                    
                    if start_location == "above"
                        [x,y] = PoleFinding.Move_xy(-dx,0,x_ll,y_ll,ARCH);
                    elseif start_location == "below"
                        [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
                    end
                    
                    %check if new pole gives line inside the arch
                    check = Toolkit.CheckThrust(ARCH,x,y);
                    
                    %if the new line is outside
                    if check == "out"
                        %calculate maximum distance from thrust line node to
                        %its closest arch geometry coordinate
                        error_new = Toolkit.max_of_diffs(ARCH,x,y);
                        
                        current_location = Toolkit.AboveOrBelow(ARCH,x,y);
                        
                        % if switched from 'above' to 'below' or vice versa
                        if current_location ~= start_location
                            
                            current_location = Toolkit.AboveOrBelow(ARCH,x,y);
                            if current_location == "above"
                                dx_back = - dx;
                            elseif current_location == "below"
                                dx_back = + dx;
                            end
                            
                            [x,y] = PoleFinding.Move_xy(dx_back,0,x_ll,y_ll,ARCH);
                            
                            %reduce step size
                            dx = dx/2;
                            error = Toolkit.max_of_diffs(ARCH,x,y);
                            dxs=dxs+1;
                            if dxs == 20 % max 20
                                %reset everything
                                dx_back = -xp;
                                dy_back = -yp;
                                [x,y] = PoleFinding.Move_xy(dx_back,dy_back,x_ll,y_ll,ARCH);
                                
                                
                                success = "no";
                                return
                            end
                            
                        else
                            error = error_new;
                        end
                    end
                    
                    
                    if check == "in"
                        xp_in = xp;
                        
                        if start_location == "above"
                            xp_out = xp + dx;
                        else
                            xp_out = xp - dx;
                        end
 
                    end
                end
                
            elseif check == "in"
                
                while check=="in"
                    [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
                    check = Toolkit.CheckThrust(ARCH,x,y);
                end
                xp_out = xp;
                xp_in = xp-dx;
                
            end
            
            
            
            
            error = Toolkit.min_of_diffs(ARCH,x,y);
            while error>tol
                %use bisection method to iterate to accurate solution
                
                xp_new = (xp_out+xp_in)/2;
                dx = xp_new - xp;
                [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
                
                check = Toolkit.CheckThrust(ARCH,x,y);
                %update either xp_in or xp_out
                if check == "out"
                    xp_out = xp;
                    error = Toolkit.max_of_diffs(ARCH,x,y);
                elseif check == "in"
                    xp_in = xp;
                    error = Toolkit.min_of_diffs(ARCH,x,y);
                end
                
            end
            
            %make sure thrust line lies within:
            xp_new = xp_in;
            dx = xp_new - xp;
            [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
            
            success = "yes";
            
        end
        
        function [x,y,success] = SW_find_target_height(ARCH,x,y,x_ll,y_ll,target_height,target_x,target_y,tol)
            %Updates position of pole to find minimum thrust line
            global xp;
            global yp;
            
            %move pole until the line is below target height
            dx = max(0.1*abs(x_ll(1)),0.05*abs(y_ll(1)));
            
            while max(y)>target_height
                [x,y] = PoleFinding.Move_xy(-dx,0,x_ll,y_ll,ARCH);
            end
            xp_below = xp;
            xp_above = xp+dx;
            
            error = 10*tol;
            iter = 1;
            while error > tol
                iter = iter+1;
                xp_new = (xp_above+xp_below)/2;
                dx = xp_new - xp;
                [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);

                if max(y)>target_height
                    xp_above = xp;
                    error = max(y)-target_height;
                elseif max(y) < target_height
                    xp_below = xp;
                    error = target_height - max(y);
                end
                
                if iter>20
                    break
                end

            end

            [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
            check = Toolkit.CheckThrust(ARCH,x,y);
            
            if check == "out"
                dx = xp_below - xp;
                [x,y] = PoleFinding.Move_xy(dx,0,x_ll,y_ll,ARCH);
            end
            
            check = Toolkit.CheckThrust(ARCH,x,y);
            
            if check == "in"
                success = "yes";
            else
                success = "no";
            end
        end
    end
end
