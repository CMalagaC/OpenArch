function block = ShiftBlock(block,loc,xh,yh)

xin = block.xin;
yin = block.yin ;
xex = block.xex ;
yex = block.yex ;

if loc == "ex"
    dx = xex(1) - xh;
    dy = yex(1) - yh;
    
    xin = xin - dx;
    yin = yin - dy;
    xex = xex - dx;
    yex = yex - dy;
elseif loc == "in"
    
    dx = xin(1) - xh;
    dy = yin(1) - yh;
    
    xin = xin - dx;
    yin = yin - dy;
    xex = xex - dx;
    yex = yex - dy;
end

block = struct;
block.xin = xin;
block.yin = yin;
block.xex = xex;
block.yex = yex;
end