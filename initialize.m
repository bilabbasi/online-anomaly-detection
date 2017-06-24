function [f,H,D] = initialize(X,varargin)
%------------------------------------------------------------------------%
%     Computes the density of the data given by X using a kernel density
%     estimator. A histogram is used (binning data). If no T is specified,
%     then the entire data set X will be used. The computational domain is
%     assumed to be the box [0,1] x [0,1].
% 
%     Parameters
%     ----------
%     X : Dyads for the data, assumed to be size (M x d)
%     T : History size. If no T is specified, T = size(X,1).
% 
%     Returns
%     -------
%     f : Density of the selected dyads.
%     H : Data being used to compute the density f.
%     D : List of dyads associated with H.
%------------------------------------------------------------------------%
grid_setup
M = size(X,1);  %Number of data points

if nargin == 1
    H = X;
else
    T = varargin{1};
    I = randperm(M,T);  %Randomly select T samples
    H = X(I,:);
    M = T;
end

D = test_dissimilarity(H);
% D = dissimilarity(H);
f = zeros(N);        %Preallocating memory for density
I = floor(D/h+1);    %Binning values to the bottom left corner of each cell
for k = 1 : size(D,1)
    f(I(k,1),I(k,2)) = f(I(k,1),I(k,2)) + 1;
end
f = f/(M*h^2); %Normalization
    
end