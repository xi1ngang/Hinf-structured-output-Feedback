% This function compute the mininorm element of the convex full of some
% poitns specificed in matrix P. This problem can be veiwed as a quadratic
% programming, and can be solved via quadprog commned in Matlab. 
% Input P: a n by m matrix where n is the dimension of the problem and m is
% the number of points
% Outpu x: a n by 1 vector which is the mininorm element in the convex hull
function x = mininorm_convexhull(P)

m = size(P,2);
H = P'*P;
lb = zeros(1,m);
Aeq = ones(1,m);
beq = 1;


options = optimoptions('quadprog','Display','none');
lambda = quadprog(H,[],[],[],Aeq,beq,lb,[],[],options);

x = P*lambda;

end