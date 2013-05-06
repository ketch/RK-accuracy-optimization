%Optimize accuracy for spectral semi-discretization of advection

nx = 50.;
dx = 1./(nx+1.);
cflnum = 0.9;
dt = cflnum * dx;

lam=spectrum('imagaxis',nx)/dx;
%xsi = (1:(nx+1))';
%lam = -1i*sin(2*pi*xsi*dx)/dx
[status,x_opt,v] = opt_acc(dt,lam,s,p);

mon_poly_coeff = [1./factorial(0:p) x_opt'];
pc=mon_poly_coeff(end:-1:1);

semilogy(imag(dt*lam),abs(polyval(pc,dt*lam)-exp(dt*lam)),'ok','linewidth',3);
hold on
pc2 = [1./factorial(s:-1:0)];
semilogy(imag(dt*lam),abs(polyval(pc2,dt*lam)-exp(dt*lam)),'or','linewidth',3);
legend('accuracy-optimal method','order-optimal method','Location','Best')
hold off

