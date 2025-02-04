function F=Objective_Surface(x,Counts,D,s)

K=zeros(size(Counts));
for kk=1:length(K)
    K(kk)=Waterfowl_Surface(D(kk,:),10.^x,s);
end

F=K(:)-Counts(:);

F=mean(F.^2);

end