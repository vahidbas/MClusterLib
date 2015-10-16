function r = mnrand_draw(p,m)

W = mnrnd(1,p,m);
for i=1:m
    r(i,1) = find(W(i,:));
end
