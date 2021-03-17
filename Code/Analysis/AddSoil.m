%% OpenArch - Function to add vertically distributed mass with an associated inertial contribution

%  Contact:
%  T. McLean, thomas.o.mclean@outlook.com
%  C. Malaga-Chuquitaype, c.malaga@imperial.com

function geom = AddSoil(geom,L,alpha_m)

global GravitationalAcceleration
global LateralAcceleration

%L is proportional multiplier of the arch's mass
%alpha_m is a factor to take into account the effective mass of the soil

%Calculate the total mass of the arch
M_arch = 0;
for i =1:length(geom)
    M_arch = M_arch + geom(i).M;
end

%Calculate the total mass of soil
M_soil = M_arch*L;

%% Distribute the mass as a udl

%extract the extrados of the arch
extrados = zeros(1,length(geom)+1);
thetas = zeros(1,length(geom)+1);

for i =1:length(geom)
    extrados(i) = geom(i).r(3);
    thetas(i) = geom(i).theta(3);
end
extrados(end) = geom(end).r(4);
thetas(end) = geom(end).theta(4);

[x,y] = Toolkit.pol2car(extrados,thetas);

%span of the arch
S = abs(x(1)-x(end));

%create vertical udl of soil mass:
w = M_soil/S;

%% Calculate loads applied to each voussoir

%vertical projection of extrados to the udl defines contribution
dx = zeros(1,length(geom));

for i = 1:length(geom)
    dx(i) = abs(x(i)-x(i+1));
end

% mass of soil applied to each voussoir
m = dx.*w;

Fv = m*GravitationalAcceleration;
Fh = alpha_m*m*LateralAcceleration;

%% Apply loads at the voussoir centroids

for i = 1:length(geom)
    geom(i).W = geom(i).W + Fv(i);
    geom(i).Fh = geom(i).Fh + Fh(i);   
end

end