function Lambda = AbateCost(mu,theta1,theta2)
%theta1 is a variable and theta2 is scale
Lambda(1,:) = theta1/1000.*mu.^theta2;