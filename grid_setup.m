N = 100;        %Number of grid points used in each coordinate
h = 1/(N-1);    %Spatial resolution
[x,y] = meshgrid(linspace(0,1,N));  %Computational domain

k1 = 10;            %Number of neighbours used for first criteria
k2 = 10;            %Number of neighbours used for second criteria