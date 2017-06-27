%Solving the equationsystem.
step=1;
intervallength=100;
intervalstart=0;
syms R T;
k=1;
acc=3;
boundryconditions=[vpa(exp(1i*intervalstart*k),acc)-R*vpa(exp(-1i*intervalstart*k),acc),T*vpa(exp(1i*intervallength*k),acc)];
derivativeconditions=[1i*k*exp(1i*intervalstart*k)+R*1i*k*exp(-1i*intervalstart*k),1i*k*T*exp(1i*intervallength*k)];
E=1;
h_bar=1;
m=1;
size=intervallength-intervalstart;
size=size/step;
interval=[intervalstart+step:step:intervallength];

v=zeros(1,size);
H=zeros(size,size);
b=zeros(size,1);
phi=zeros(size,1);
eigenvectors=zeros(size);
eigenvalues=zeros(size);

%Input boundry derivative data into b vector
b=[boundryconditions(1);zeros(size-2,1);boundryconditions(2);];
%b(1)=boundryconditions(1);
%b(size)=boundryconditions(2);

%Creating the potential
%v(5/step:15/step)=0.1*sin(0:pi/(5/step):2*pi);
v(40/step:60/step)=1;

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
          H(i,j)=-1/(step^2)*h_bar/(2*m);
       end
       if j==i
       H(i,j)= 2/(step^2)*h_bar/(2*m)+v(j)-E;
       end
   end
end
%H(size+1,1)=1/step;
%H(size+1,2)=-1/step;
%H(size+2,size-1)=1/step;
%H(size+2,size)=-1/step;

%answer=linsolve(H,b);
%Get a numerical solution to the problem H*phi =b
%answer=inv(H)*b;
answer=H\b;
%Now find R,T from this
temp1=solve(answer(1)/step-answer(2)/step == derivativeconditions(1), answer(size)/step-answer(size-1)/step == derivativeconditions(2),R,T);
temp2=[temp1.R temp1.T];
R=temp2(1);
T=temp2(2);
answer=subs(answer);
answer=abs(answer);
subplot(2,1,1), plot(interval',answer)
subplot(2,1,2), plot(interval',v)

%{
Problems:
*R not 0 and T not 1 when there is no potential present!
*Very slow solver when step is not 1. 


%}


