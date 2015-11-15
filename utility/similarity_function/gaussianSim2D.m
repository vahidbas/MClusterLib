function s = gaussianSim2D(x,y)

B = eye(2); % covarience matrix

s = 1/(2*pi*sqrt(det(B)))*exp(-0.5*(x-y)/B*(x-y)');