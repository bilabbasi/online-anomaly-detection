function dyads = norm_pt(dyads,D)
    %Returns dyads normalized wrt D (component-wise).
    Dmax = max(D);
    Dmin = min(D);
    if size(dyads,1) == 1
        dyads = max(min((dyads-Dmin)./(Dmax-Dmin),1),0);
    else
        for i = 1 : length(dyads)
           dyads(i,:) = max(min((dyads(i,:)-Dmin)./(Dmax-Dmin),1),0);
        end
    end
end