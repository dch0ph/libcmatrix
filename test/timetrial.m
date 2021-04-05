sys=[1/2,1/2,1/2];

H01=10*spindipolar(sys,1,2)+15*spindipolar(sys,2,3);
H02=100*spindipolar(sys,1,3);

UNIT=100000;
UNITST=10000;
n=8;

CADDS=UNIT*n*n;

b=2+3*i;
c=4+3*i;

disp('Complex number multiplication');

t=cputime;
for k=1:CADDS,
  a=b*c;
end
used=cputime-t;
disp(used);
disp(a);


disp('Matrix addition');

t=cputime;
for k=1:UNIT
  H=H01+H02;
end
used=cputime-t;
disp(used);
disp(H);

disp('Multiplication');
t=cputime;
for k=1:UNITST
  H=H01*H02;
end
used=cputime-t;
disp(used);
disp(H);

U=expm(-i*1e-3*H02);
disp('Similarity transform');
disp(U);

t=cputime;
for k=1:UNITST
  H=U'*H01*U;
end
used=cputime-t;
disp(used);
disp(H);
