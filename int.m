syms S I R x e m A B a b n i d g r k
f1 = @(S,I,R,x)e*m*x.*(S+I+R)./(A*x.+B)-a*b*S.*I./(S+I+R).-n*S.-d*S.;
f2 = a*b*S*I/(S+I+R)-i*I-g*I;
f3 = d*S+g*I-n*R;
f4 = r*x*(1-x/k)-m*x*(S+I+R)/(A*x+b);

up = {10,10,10,100};
low = {0,0,0,0};
f = nIntegrate(f1,low,up)