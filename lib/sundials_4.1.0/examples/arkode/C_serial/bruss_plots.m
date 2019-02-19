% ----------------------------------------------------------------
% Programmer(s): Daniel R. Reynolds @ SMU
% ----------------------------------------------------------------
% SUNDIALS Copyright Start
% Copyright (c) 2002-2019, Lawrence Livermore National Security
% and Southern Methodist University.
% All rights reserved.
%
% See the top-level LICENSE and NOTICE files for details.
%
% SPDX-License-Identifier: BSD-3-Clause
% SUNDIALS Copyright End
% ----------------------------------------------------------------
% Matlab script to load/plot 1D brusselator results

clear

% load output files
x = load('bruss_mesh.txt');
u = load('bruss_u.txt');
v = load('bruss_v.txt');
w = load('bruss_w.txt');

% check for compatible sizes
Nx = length(x);
[Nt_u,Nx_u] = size(u);
[Nt_v,Nx_v] = size(v);
[Nt_w,Nx_w] = size(w);
if ((Nx ~= Nx_u) | (Nx ~= Nx_v) | (Nx ~= Nx_w))
   fprintf('Error: incompatible spatial dimensions!\n');
   fprintf('  Nx = %i,  Nx_u = %i,  Nx_v = %i,  Nx_w = %i\n',Nx,Nx_u,Nx_v,Nx_w);
   error('Result files are incompatible, re-run simulation')
end
Nt = Nt_u;
if ((Nt ~= Nt_v) | (Nt ~= Nt_w))
   fprintf('Error: incompatible temporal dimensions!\n');
   fprintf('  Nt_u = %i,  Nt_v = %i,  Nt_w = %i\n',Nt_u,Nt_v,Nt_w);
   error('Result files are incompatible, re-run simulation')
end

% get bounding box size
ymax = 1.1*max([max(max(u)), max(max(v)), max(max(w))]);
ymin = 0.9*min([min(min(u)), min(min(v)), min(min(w))]);
xmax = max(x);
xmin = min(x);

% plot time series for center of spatial domain
figure()
t = linspace(1,Nt,Nt);
ix = floor(Nx/2);
plot(t,u(:,ix),'b-',t,v(:,ix),'r--',t,w(:,ix),'k-.','LineWidth',1)
legend('u','v','w')
axis([1, Nt, ymin, ymax]);
title('Time series of solutions at domain center')


% loop over output times, plotting results and pausing between each
figure()
for it = 1:Nt
   uvec = u(it,:)';
   vvec = v(it,:)';
   wvec = w(it,:)';
   plot(x,uvec,'b-',x,vvec,'r--',x,wvec,'k-.','LineWidth',1)
   legend('u','v','w')
   axis([xmin, xmax, ymin, ymax]);
   title(sprintf('Solutions, time step %i',it-1))
   pause
end
