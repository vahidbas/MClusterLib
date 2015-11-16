
x = rand(1,1000)*20-10;
y = rand(1,1000)*20-10;

for i=1:length(x)
    X{i} = [x(i);y(i)];
end

alg = makeClustring('itm',3,'Plot','on');
result = alg.cluster(X);
