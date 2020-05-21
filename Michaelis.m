syms v Ubv Ki I S Km

Vm=10.0;
Km=21.0;

Ki= [3722.727 2863.636364 13077.27273 763.6363636 ];    %Ki for DUB
Kig= [10309.09091 15654.54545 84859.09091 1145.454545]; %Ki for deISGylation

S=1.0;


I=linspace(0.1,10000,100);
for j=1:4
for i=1:100
    

v= (Vm*S/((1+((I(i))/Ki(j))*Km)+(S)));
v2= (Vm*S)/((1+((I(i))/Kig(j))*Km)+(S));


Soln(i)=v;

Soln2(i)=v2;
end
hold on
figure (1)
plot( log10(I), Soln);

figure (2)
plot( log10(I), Soln2);




end 




