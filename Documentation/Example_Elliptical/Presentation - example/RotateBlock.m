function block = RotateBlock(block,phi,loc,xh,yh)
%phi angle of rotation
%loc - location of lower hinge (extrados or intrados)
%xh and yh give the location of the hinge

xin = block.xin;
yin = block.yin ;
xex = block.xex ;
yex = block.yex ;

if loc == "in"
    dx = xin(1);
    
    xin = xin - dx;
    xex = xex - dx;
    
    T = [cos(phi) -sin(phi);
        sin(phi) cos(phi)];
    
    for i = 1:length(xin)
        X = T*[xin(i);yin(i)];
        xin(i) = X(1);
        yin(i) = X(2);
        
        X = T*[xex(i);yex(i)];
        xex(i) = X(1);
        yex(i) = X(2);
    end
    
    
    xin = xin + dx;
    xex = xex + dx;
    
    
    block = struct;
    block.xin = xin;
    block.yin = yin;
    block.xex = xex;
    block.yex = yex;
    
    
elseif loc == "ex"
    
    dx = xex(1);
    
    xin = xin - dx;
    xex = xex - dx;
    
    T = [cos(phi) -sin(phi);
        sin(phi) cos(phi)];
    
    for i = 1:length(xin)
        X = T*[xin(i);yin(i)];
        xin(i) = X(1);
        yin(i) = X(2);
        
        X = T*[xex(i);yex(i)];
        xex(i) = X(1);
        yex(i) = X(2);
    end
    
    
    xin = xin + dx;
    xex = xex + dx;
    
    
    block = struct;
    block.xin = xin;
    block.yin = yin;
    block.xex = xex;
    block.yex = yex;
    

end
end