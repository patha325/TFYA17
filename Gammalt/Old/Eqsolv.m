%Solving the equationsystem.
step=1;
intervallength=20;
intervalstart=0;
boundryconditions=[1,1];
energy=1;
h_bar=1;
m=1;
size=intervallength-intervalstart;
size=size/step;
interval=[intervalstart+step:step:intervallength];

%interval=intervalstart:step:intervalstart+intervallength;

v=zeros(1,size);
H=zeros(size);
b=zeros(size,1);
phi=zeros(size,1);
eigenvectors=zeros(size);
eigenvalues=zeros(size);

%Input boundry data into phi vector
%phi(1)=boundryconditions(1);
%phi(size)=boundryconditions(2);

%Creating the potential
v(5/step:15/step)=1;


%Get the eq. H*phi =b
%Input boundry into H
H(1,1)=2/(step^2)*h_bar/(2*m)+v(1);
H(1,2)=-1/(step^2)*h_bar/(2*m);
%H(1,3)=-1/(step^2)*h_bar/(2*m);
H(size,size)=2/(step^2)*h_bar/(2*m)+v(size);
H(size,size-1)=-1/(step^2)*h_bar/(2*m);
%H(size,size-2)=-1/(step^2)*h_bar/(2*m);

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

%Create b
for i=1:size,
    b(i)=energy;
end

%Eigenvalues & eigenvectors to H.
[eigenvectors,eigenvalues]=eig(H);
answer=eigenvalues*eigenvectors;
%Normalize???
subplot(2,1,1), plot(interval',answer(:,1:5))
subplot(2,1,2), plot(interval',v)
%plot(interval',eigenvectors(:,1:5))

