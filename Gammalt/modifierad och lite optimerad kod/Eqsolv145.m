function conductance =Eqsolv14(E,step)

%Solving the equationsystem.
%step=1;
intervallength=10;
intervalstart=1;
syms R T;
%E=1;
h_bar=1;
m=1;
k=sqrt(2*E*m/h_bar^2);
acc=5;
boundryconditions=[vpa(exp(1i*intervalstart*k),acc)-R*vpa(exp(-1i*intervalstart*k),acc),T*vpa(exp(1i*intervallength*k),acc)];
derivativeconditions=[1i*k*exp(1i*intervalstart*k)+R*1i*k*exp(-1i*intervalstart*k),1i*k*T*exp(1i*intervallength*k)];
size=intervallength-intervalstart;
size=size/step;
interval=[intervalstart+step:step:intervallength];

v=zeros(1,size);
H=zeros(size,size);
b=zeros(size,1);
phi=zeros(size,1);

%Input boundry derivative data into b vector
b=[boundryconditions(1);zeros(size-2,1);boundryconditions(2);];

beta=1/2;
Vsd=0.3;
Vg=0.3;
s=5;
x1=4;
x2=6;
x=intervalstart:step:intervallength-1*step;
%Creating the potential
%v(5/step:15/step)=0.1*sin(0:pi/(5/step):2*pi);
%v(4/step:6/step)=1;
v=(beta-1/2)*Vsd+Vg*tanh(s*(x-x1))-(Vsd/2+Vg)*tanh(s*(x-x2));
%v=0.5*tanh((intervalstart/step:intervallength/step-1)-4/step)-0.3*tanh((intervalstart/step:intervallength/step-1)-6/step);


%Get the eq. H*phi =b
%Input boundry into H
%H(1,1)=2/(step^2)*h_bar/(2*m)+v(1);
%H(1,2)=-1/(step^2)*h_bar/(2*m);
%H(size,size)=2/(step^2)*h_bar/(2*m)+v(size);
%H(size,size-1)=-1/(step^2)*h_bar/(2*m);
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
%H(size+1,1)=1/step;
%H(size+1,2)=-1/step;
%H(size+2,size-1)=1/step;
%H(size+2,size)=-1/step;

answer=linsolve(H,b);
%Get a numerical solution to the problem H*phi =b
%answer=inv(H)*b;
%answer=H\b;
%Now find R,T from this

temp1=solve(answer(2)/step-answer(1)/step == derivativeconditions(1), answer(size)/step-answer(size-1)/step == derivativeconditions(2),R,T);
temp2=temp1;%[temp1.R temp1.T];
R=temp2.R(1);
T=temp2.T(1);
R=abs(R)^2/(abs(R)^2+abs(T)^2);
T=abs(T)^2/(abs(R)^2+abs(T)^2)
conductance = T/(R+T);
answer=subs(answer);
answer=abs(answer);
answer2=answer.^2;
answer2=answer2/norm(answer);
%conductance=abs(T)/abs(R);

subplot(2,1,1), plot(interval',answer2)
subplot(2,1,2), plot(interval',v)

end


%{
Problems:
*R not 0 and T not 1 when there is no potential present!
*Very slow solver when step is not 1. 


%}


