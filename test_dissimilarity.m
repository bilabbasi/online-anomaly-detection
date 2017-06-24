function D = test_dissimilarity(H,varargin)
%------------------------------------------------------------------------%
%     Dissimilarity example used for testing. Dissimilarity metric here is
%     just component wise difference:
%                  d(a,b) = [abs(a1-b1), abs(a2-b2)]
%     
%     Parameters
%     ----------
%     H : Old data samples used. Assumed to be of the size Mx2.
%
%     Returns
%     -------
%     D : dyads                                           
%------------------------------------------------------------------------% 
if nargin == 1
    T = size(H,1);
    num_dyads = nchoosek(T,2);
    D = zeros(num_dyads,2);
    k = 1;
    for i = 1 : T
        for j = i+1 : T
            D(k,:) = [abs(H(i,1)-H(j,1)) abs(H(i,2)-H(j,2))];
            k = k + 1;
        end
    end
else %If x point is also given in addition to H then only compute dyads for x
    x = cell2mat(varargin);
    D = zeros(size(H,1),2);
    for k = 1 : size(H,1)
       D(k,:)  = [abs(x(1) - H(k,1)) abs(x(2) - H(k,2))];
    end
end

end