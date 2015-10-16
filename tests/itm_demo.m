addpath('./../ITM') % add ITM class directory

x = rand(1,1000)*20-10;
y = rand(1,1000)*20-10;

for s = 1:5;

% x = Xm{s}+randn(size(Xm{s}));
% y = Ym{s}+randn(size(Xm{s}));

% x = Xm{s};
% y = Ym{s};



b =1;
if s == 1
 net = ITM(3);
 b = 3;
end

for t = b:length(x) 
    
    step(net, [x(t);y(t)],eye(2));
    if net.num > 2
    plot(net);
    hold on
    plot(x(t),y(t),'ob');
    hold off
    end
    pause(0.1)
end
end

for k = 1:net.num
    for m = 1:net.num
        dmk(m,k) = sqrt(norm(net.nodes(k).w-net.nodes(m).w))/2;
    end
end
