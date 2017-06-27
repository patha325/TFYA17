function dhdt = tankfill(t,h)
% RHS function for tank-fill problem
A = 10 + 4*sin(t);
% alpha(t)
H = 2*sqrt(h);
% beta*sqrt(h)
dhdt = A - H;
% eof- tankfill.m