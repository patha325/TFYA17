% In here is the solution to the schrödinger eq given in TFYA17fkn. 
x0=[1 1];
tspan=[0,20];
[t,x]=ode45('TFYA17fkn',tspan,x0);
plot(t,x(:,1))
