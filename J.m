syms S I R x e m A B a b n i d g r k L C0 C1 C2 C3 C4 C5
f1 = e*m*x*(S+I+R)/(A*x+B)-a*b*S*I/(S+I+R)-n*S-d*S;
f2 = a*b*S*I/(S+I+R)-i*I-g*I;
f3 = d*S+g*I-n*R;
f4 = r*x*(1-x/k)-m*x*(S+I+R)/(A*x+b);

f1S = diff(f1,S);
f2S = diff(f2,S);
f3S = diff(f3,S);
f4S = diff(f4,S);
f1I = diff(f1,I);
f2I = diff(f2,I);
f3I = diff(f3,I);
f4I = diff(f4,I);
f1R = diff(f1,R);
f2R = diff(f2,R);
f3R = diff(f3,R);
f4R = diff(f4,R);
f1x = diff(f1,x);
f2x = diff(f2,x);
f3x = diff(f3,x);
f4x = diff(f4,x);

J0 = [f1S f1R f1I f1x;
     f2S f2R f2I f2x;
     f3S f3R f3I f3x;
     f4S f4R f4I f4x
    ];
I1 = eye(4);

% J1 = J0;
% for i = 1:4
%     for j = 1:4
%         J1(i,j) = subs(J1(i,j),{S,I,R,x},{C0*(C1+1)*C2*(A*C5+B)*(1-C5/k),C1*C2*(A*C5+B)*(1-C5/k),C2*(A*C5+B)*(1-C5/k),C5});
%     end
% end
% LJ1 = J1 - L*I1;
% D = det(LJ1);
% eq = D==0;
% s = solve(eq,L)
% solL = roots(s)

J2 = J0;
for i = 1:4
    for j = 1:4
        J2(i,j) = subs(J2(i,j),{S,I,R,x},{e*r*((B*C2)/(m*e*(C1+C2)-A*C1*C2)-(B^2*C1*C2^2)/([m*e*(C1+C2)-A*C1*C2]^2*k)),0,e*r*((B*C1)/(m*e*(C1+C2)-A*C1*C2)-(B^2*C1^2*C2)/([m*e*(C1+C2)-A*C1*C2]^2*k)),(B*C1*C2)/(m*e*(C1+C2)-A*C1*C2)});
    end
end
LJ2 = J2 -L*I1;
D = det(LJ1);
eq = D==0;
s = solve(eq,L)
solL = roots(s)