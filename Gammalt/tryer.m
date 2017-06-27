function tryer(steps,res,res2)
c=zeros(1,steps);
for i=0:steps/res
    c(1,i+1)=Eqsolv14(i*res,res2);
end
plot(0:steps/res,c)
end