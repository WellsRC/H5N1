function K=Waterfowl_Surface(D,Scale_K,s)
K = sum(Scale_K(:).*normpdf(D(:),0,s)./normpdf(0,0,s));
end