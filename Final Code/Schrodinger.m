function conductance =Schrodinger(E,step,s,color)
%Schrodinger(E,steps,s,color) plots the solution to the Schrodinger equation, the input
%potential and also calculates the conductance of the system. The program
%also requires Matlab 2013a or newer.
%Input the following:
%Energy that the incoming particle has in CGS units.
%The size of steps to divide the interval in, default interval is 1:10
%steps is recomended to be 0.1-0.01 any larger will increase the calculation time without improving the calculations significantly.
%The value of s in the potential, since the default potential is 
%v=(beta-1/2)*Vsd+Vg*tanh(s*((intervalstart/step:intervallength/step-1)-x1/step))-(Vsd/2+Vg)*tanh(s*((intervalstart/step:intervallength/step-1)-x2/step));
%Vsd and Vg are set internaly to 0.3, also beta is set to 1/2.
%The colour of the plot, given as e.g. 'k' for black, 'r' for red.
%
hold on
%Solving the equationsystem.
intervallength=10;
intervalstart=1;
syms R T;
h_bar=1;
m=1;
e=1; %Electron charge
beta=1/2;
Vsd=0.3;
Vg=0.3;
x1=4;
x2=6;

k=sqrt(2*m/h_bar^2 * (E-beta*e*Vsd));
q=sqrt(2*m/h_bar^2 * (E+(1-beta)*e*Vsd));


acc=5;
boundryconditions=[vpa(exp(1i*intervalstart*k),acc)+R*vpa(exp(-1i*intervalstart*k),acc),T*vpa(exp(1i*intervallength*q),acc)];
derivativeconditions=[1i*k*exp(1i*intervalstart*k)-R*1i*k*exp(-1i*intervalstart*k),1i*q*T*exp(1i*intervallength*q)];
size=intervallength-intervalstart;
size=size/step;
interval=[intervalstart+step:step:intervallength];
v=zeros(1,size);
H=zeros(size,size);
b=zeros(size,1);
phi=zeros(size,1);
%Input boundry derivative data into b vector
b=[boundryconditions(1);zeros(size-2,1);boundryconditions(2);];


x=intervalstart:step:intervallength-1*step;
%Creating the potential
%v(5/step:15/step)=0.1*sin(0:pi/(5/step):2*pi); Sinus
%v(4/step:6/step)=1; Simple step.
%v=0.5*tanh((intervalstart/step:intervallength/step-1)-4/step)-0.3*tanh((intervalstart/step:intervallength/step-1)-6/step);
v=(beta*e-1/2)*Vsd+Vg*tanh(s*((intervalstart/step:intervallength/step-1)-x1/step))-(Vsd/2+Vg)*tanh(s*((intervalstart/step:intervallength/step-1)-x2/step));

H(1,1)=1;
H(size,size)=1;

%Calculate the H matrix. 
for i=2:size-1,
   for j=1:size, 
       if j==i-1 || j==i+1,
          H(i,j)=-h_bar^2/(2*m*step^2);
       end
       if j==i
       H(i,j)= h_bar^2/(m*step^2)+v(j)-E;
       end
   end
end
%Get a numerical solution to the problem H*phi =b
answer=linsolve(H,b);

%Now find R,T from this
temp1=solve(answer(2)/step-answer(1)/step == derivativeconditions(1), answer(size)/step-answer(size-1)/step == derivativeconditions(2),R,T);
temp2=temp1;%[temp1.R temp1.T];
R=temp2.R(1);
T=temp2.T(1);
R=abs(R)^2/(abs(R)^2+abs(T)^2);
T=abs(T)^2/(abs(R)^2+abs(T)^2);
conductance = T/(R+T);
answer=subs(answer);
answer=abs(answer);
answer2=answer.^2;
answer2=answer2/norm(answer);
%conductance=abs(T)/abs(R);

subplot(2,1,1)
hold on
plot(interval',answer2,color)
ylabel('|\Psi (x)|^2')
xlabel('X')
subplot(2,1,2)
hold on
plot(interval',v,color)
ylabel('V(X)')
xlabel('X')

end