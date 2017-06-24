function [anomaly_score,anomaly_class] = anomaly_score(f,H,D,p,threshold)
%------------------------------------------------------------------------%
%     Uses PDE-based sorting to compute the anomaly scores and consequently
%     determine whether or not the points in x are anomalies. Anomalous
%     points are further classified based on the type of anomaly.
% 
%     Parameters
%     ----------
%     f : Estimated density of D 
%     H : Old data samples used
%     D : Old dyads used to construct f. Needed to normalize the new dyads.
%     p : New data
%     threshold : Threshold for anomaly scoring. If anomaly score is bigger
%                 than the treshold than the data is an anomaly.
%
%     Returns
%     -------
%     anomaly_score : Binary classification whether point is anomaly or
%                     not. 
%                     1 = anomaly and 0 = normal
%     anomaly_class : Binary classification of which criteria is violated.
%                     1 = first, 2 = second, 3 = inconclusive
%                     0 = not an anomaly
%------------------------------------------------------------------------%

grid_setup

U = sortPDE(f+h^2);       %Compute scoring function
W = classPDE(U,f+h^2);    %Compute classifying function
anomaly_score = zeros(size(p,1),1);
anomaly_class = zeros(size(p,1),1);

for i = 1 : size(p,1)
    
    pdyads = test_dissimilarity(H,p(i,:));
    pdyads = norm_pt(pdyads,D);
    [~,k1index] = sort(pdyads(:,1));
    [~,k2index] = sort(pdyads(:,2));
    nbrs = [pdyads(k1index(1:k1),:); 
            pdyads(k2index(1:k2),:)];
    pscore = mean(interp2(x,y,U,nbrs(:,1),nbrs(:,2)));
    if pscore > threshold
        anomaly_score(i) = 1;
        pclass = mean(interp2(x,y,W,nbrs(:,1),nbrs(:,2)));
        if pclass > 0.55
            anomaly_class(i) = 1;
        elseif pclass < 0.45
            anomaly_class(i) = 2;
        else
            anomaly_class(i) = 3;
        end
    end

end
end

function U = sortPDE(F)
%% Description
% Returns the numerical solution to U_xU_y=F by constructing a solution to
% the scaled PDE W_xW_y=WF, where W = U^2/4. The construction uses a first
% order backwards difference near the boundary, and uses a second order
% backwards difference in the rest of the interior. 
%% Variables
% F - Density of sample used in sorting.
% U0 - Initial Gauss-Seidel guess.
% N - Grid size, i.e., solution is given on an NxN grid.
% h - Step size.
% W - Solution of scaled PDE.
%% Construction of numerical solution
N = size(F,1);
h = 1/(N-1);
W = zeros(N);

for i = 2 : N
    for j = 2 : N
        W(i,j) = 0.5*(W(i-1,j) + W(i,j-1) + h^2*F(i,j)) + 0.5*sqrt(...
            (W(i,j-1)-W(i-1,j))^2+2*h^2*F(i,j)*(W(i-1,j)+W(i,j-1))+...
            h^4*F(i,j)^2);
    end
end
U = 2*W.^(1/2); %Scales back to original solution
end

function [W,V] = classPDE(U,F)
%% Description
% Returns the solution to the PDE given by UyWx-UxWy=FW/V.

V = classPDEaux(U,F);
N = size(U,1);
h = 1/(N-1); %Resolution

% Boundary data of 1 on right and bottom boundaries. Downwind sweep from
% bottom right to top left.
W = ones(N);
for i = 2 : 1 : N
    for j = N-1 : -1 : 1
        Ux = (U(i,j+1)-U(i,j))/h;
        Uy = (U(i,j)-U(i-1,j))/h;
        W(i,j) = (Uy*W(i,j+1)+Ux*W(i-1,j))/(Ux+Uy+h*F(i,j)/V(i,j));
    end
end
W(end,:) = 0;
W(:,1) = 0;
end

function V = classPDEaux(U,F)
%% Description
% Returns the solution to the PDE given by UyVx-UxVy=F, with zero boundary
% data on the top and left boundaries. This solution is then used to
% compute the normalized version used for the actual anomaly
% classification.
%% Variables
% F - density
% U - approximate ranking function
% V - solution to PDE (V(x,y) returns number of points above (x,y) in the
% front to which (x,y) belongs.
N = size(U,1);
h = 1/(N-1);
V = zeros(N);
for i = N-1 : -1 : 2
    for j = 2 : 1 : N
        Ux = (U(i,j)-U(i,j-1))/h;
        Uy = (U(i,j)-U(i-1,j))/h;
        
        if Ux + Uy > 0
            V(i,j) = (Uy*V(i,j-1) + Ux*V(i+1,j) + h*F(i,j))/(Ux+Uy);
        else
            V(i,j) = (V(i,j-1) + V(i+1,j))/2;
        end
    end
end

end
