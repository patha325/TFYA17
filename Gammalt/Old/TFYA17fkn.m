% In here is the schrödinger eq. 

function phi=schrodinger(t,x);
h_bar=1;
m=1;
v=t;
phi=zeros(size(x));
phi(1)=x(2);
phi(2)=v*x(1)*2*m/(h_bar^2);






