function [h,k,error] = spectral_advection(m,pc)
%
% Solve u_t + au_x = 0  on [ax,bx] with periodic boundary conditions,
% using spectral collocation in space and an ERK (polynomial) in time.
%
% Returns k, h, and the max-norm of the error.
% This routine can be embedded in a loop on m to test the accuracy,
% perhaps with calls to error_table and/or error_loglog.
%

global a
a = 1;           % advection velocity

clf              % clear graphics

ax = 0;
bx = 1;
tfinal = 1;                % final time

h = (bx-ax)/(m);         % h = delta x
nu = 1.;                % Courant number

k = nu*h/a;
x = linspace(ax,bx,m+1)';  % note x(1)=0 and x(m+2)=1
x = x(1:end-1);
                           % With periodic BC's there are m+1 unknowns u(2:m+2)
I = 2:(m+2);   % indices of unknowns

nsteps = round(tfinal / k);    % number of time steps
nplot = 20;       % plot solution every nplot time steps
                  % (set nplot=2 to plot every 2 time steps, etc.)
%nplot = nsteps;  % only plot at final time

L = bx-ax;
xsi = (2*pi/L)*[0:(m/2-1) (-m/2):-1]'; % Wavenumber vector in either x or y

if abs(k*nsteps - tfinal) > 1e-5
   % The last step won't go exactly to tfinal.
   disp(' ')
   disp(sprintf('WARNING *** k does not divide tfinal, k = %9.5e',k))
   disp(' ')
   end

% initial conditions:
tn = 0;
u0 = eta(x);
u = u0;

%%%%%%%%%%%%%%%%%%%%%%
% Optimized method
%PA = polyvalm(pc,-1i*xsi*a*k);
PA = polyval(pc,-1i*xsi*a*k);
%%%%%%%%%%%%%%%%%%%%%%

% initial data on fine grid for plotting:
xfine = linspace(ax,bx,1001);
ufine = utrue(xfine,0);

% plot initial data:
plot(x,u0,'b.-', xfine,ufine,'r')
axis([0 1 -.2 1.2])
legend('computed','true')
title('Initial data at time = 0')

input('Hit <return> to continue  ');

% main time-stepping loop:

for n = 1:nsteps
     tnp = tn + k;   % = t_{n+1}

    uhat = fft(u);
    %uhat = exp(-1i*xsi*a*k).*uhat;
    uhat = PA.*uhat;
    u = ifft(uhat);



     % plot results at desired times:
     if mod(n,nplot)==0 | n==nsteps
        uint = u(1:m);  % points on the interval (drop ghost cell on right)
        ufine = utrue(xfine,tnp);
        plot(x,uint,'b.-', xfine,ufine,'r')
        axis([0 1 -.2 1.2])
        title(sprintf('t = %9.5e  after %4i time steps with %5i grid points',...
                       tnp,n,m+1))
        error = max(abs(uint-utrue(x,tnp)));
        disp(sprintf('at time t = %9.5e  max error =  %9.5e',tnp,error))
        if n<nsteps, input('Hit <return> to continue  '); end;
        end

     tn = tnp;   % for next time step
     end
%--------------------------------------------------------

function utrue = utrue(x,t)
% true solution for comparison
global a

% For periodic BC's, we need the periodic extension of eta(x).
% Map x-a*t back to unit interval:

xat = rem(x - a*t, 1);
ineg = find(xat<0);
xat(ineg) = xat(ineg) + 1;
utrue = eta(xat);
return


%--------------------------------------------------------

function eta = eta(x)
% initial data

%beta = 600;
%eta = exp(-beta*(x - 0.5).^2);
eta = (x<0.5).*(x>0.2);
return
