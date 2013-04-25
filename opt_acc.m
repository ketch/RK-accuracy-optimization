function [status,x_opt,v] = opt_acc(h,lam,s,p)
%function [status,x_opt,v] = opt_acc(h,lam,s,p)
% Optimize accuracy (bisection subproblem)
scaling=factorial(0:s);
precision='best';

% Compute z_i^n for n=[0:s], z_i=h*lam_i
for i=1:length(lam)
     c(i,:) = (h*lam(i)).^[0:s]./scaling;
end
tall_tree_coeffs = scaling(1:p+1)./factorial(0:p);

fixedvec = c(:,1:p+1)*tall_tree_coeffs';

% Start cvx for feasible solution search 
cvx_begin
  cvx_precision(precision)
  cvx_solver('sdpt3')

  variable x_opt(s-p) 

  R=(fixedvec+c(:,p+2:end)*x_opt);
  abs(R)<=1.;
  minimize max(abs(R-exp(h*lam)));

cvx_end

status=cvx_status;
v = cvx_optval;
cvx_slvtol;

R=(fixedvec+c(:,p+2:end)*x_opt);
x_opt
scaling
length(scaling)
x_opt=x_opt./scaling(p+2:end)';
