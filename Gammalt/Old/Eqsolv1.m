%Solving the equationsystem.
step=1;
intervallength=20;
intervalstart=0;
R=sym('R');
T=sym('T');
k=1;
%boundryconditions=[R,T];
boundryconditions=[exp(1i*intervalstart*k)-R*exp(-1i*intervalstart*k),T*exp(1i*intervallength*k)];
energy=1;
h_bar=1;
m=1;
size=intervallength-intervalstart;
size=size/step;
interval=[intervalstart+step:step:intervallength];

v=zeros(1,size);
H=zeros(size);
b=zeros(size,1);
phi=zeros(size,1);
eigenvectors=zeros(size);
eigenvalues=zeros(size);

%Input boundry derivative data into b vector
b=[boundryconditions(1);zeros(size-2,1);boundryconditions(2)];
%b(1)=boundryconditions(1);
%b(size)=boundryconditions(2);

%Creating the potential
v(5/step:15/step)=100*sin(0:pi/(5/step):2*pi);


%Get the eq. H*phi =b
%Input boundry into H
H(1,1)=2/(step^2)*h_bar/(2*m)+v(1);
%H(1,2)=-1/(step^2)*h_bar/(2*m);
H(size,size)=2/(step^2)*h_bar/(2*m)+v(size);
%H(size,size-1)=-1/(step^2)*h_bar/(2*m);


%Calculate the H matrix. 
for i=2:size-1,
   for j=1:size, 
       if j==i-1 || j==i+1,
          H(i,j)=-1/(step^2)*h_bar/(2*m);
       end
       if j==i
       H(i,j)= 2/(step^2)*h_bar/(2*m)+v(j);
       end
   end
end
precition=3;
answer=inv(H)*b;
answer=vpa(answer,precition);


%subplot(2,1,1), plot(interval',answer(:,1:5))
%subplot(2,1,2), plot(interval',v)
%plot(interval',eigenvectors(:,1:5))

