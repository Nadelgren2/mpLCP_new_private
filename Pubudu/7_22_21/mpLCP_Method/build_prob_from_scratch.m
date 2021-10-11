pkg load symbolic;
syms lamda;

Q1 = [6 -6 6; -6 14 -10; 6 -10 8];
Q2 = [2 -4 2; -4 16 -2; 2 -2 4];
p1 = [0; -3; 0];
p2 = [-1; -1; 1];

A = [1 0 0;
     0 1 0;
     0 -1 0];
     
b = [3; 2; -2];

M = sym(zeros(6,6),'f');
M(1:3,4:6) = -1*A; 
M(4:6,1:3) = A';
M(4:6,4:6) = lambda*Q1 + (1-lambda)*Q2;

P = lambda*p1 + (1-lambda)*p2;

q = sym(zeros(6,1),'f');
q(1:3,1) = b;
q(4:6,1) = P;

G = sym(zeros(6,13),'f');
G(1:6,1:6) = eye(6);
G(1:6,7:12) = -1*M;
G(1:6,13) = q;
