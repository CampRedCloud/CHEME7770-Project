syms epsm eps c EC n

n=1.5;
EC= [700 900 4500 200];
epsm=0.8;

C=linspace(0.001,1000,100);
for j=1:4
for i=1:100
   
eps= (epsm*(C(i)^1.5))/(EC(j)+(C(i)^1.5));
Soln(i)=eps;

end
hold on
plot(C,Soln)

end 
