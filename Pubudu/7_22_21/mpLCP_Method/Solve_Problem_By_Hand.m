% Load symbolic package

pkg load symbolic;

%Original Data

h = 6;
k = 1;

M_data = [1 4 0 -1.00;
    2 5  0 -1;
    3 5 0 1;
    
4 1 0 1.00;
4 4 0 1.00;
4 4 1 2.00;
4 5 0 -2.00;
4 5 1 -1.00;
4 6 0 1.00;
4 6 1 2.00;

5 2 0 1;
5 3 0 -1;
5 4 0 -2.00;
5 4 1 -1.00;
5 5 0 8.00;
5 5 1 -1.00;
5 6 0 -1.00;
5 6 1 -4.00;

6 4 0 1.00;
6 4 1 2.00;
6 5 0 -1.00;
6 5 1 -4.00;
6 6 0 2.00;
6 6 1 2.00;
];

q_data = [1 0 3;
    2 0 2;
    3 0 -2;
    
4 0 -1.00;
4 1 1.00;
5 0 -1;
5 1 -2.00;
6 0 1.00;
6 1 -1.00;
];

S_theta = [1 1 -1;
2 1 1;
];

S_theta_RHS = [0;
1;
];


% Build Problem Formulation

x = sym('x',[1 k+2]);
M = sym(zeros(h,h));
q = sym(zeros(h,1));
M_data_rows = size(M_data,1);
q_data_rows = size(q_data,1);
c_data_rows = size(S_theta,1);
S_theta_rows = S_theta(end,1);
c = sym(zeros(1,S_theta_rows));
indices = zeros(1,h);

regions = [];

for i = 1:M_data_rows
    if(M_data(i,3) == 0) M(M_data(i,1),M_data(i,2)) = M(M_data(i,1),M_data(i,2)) + M_data(i,4);
    else M(M_data(i,1),M_data(i,2)) = M(M_data(i,1),M_data(i,2)) + x(M_data(i,3))*M_data(i,4);
    end
end

for i = 1:q_data_rows
%     q_data(i,:)
    if(q_data(i,2) == 0) 
        if(q_data(i,3) <= 0) q(q_data(i,1)) = q(q_data(i,1)) + q_data(i,3) + x(k+1)*(abs(q_data(i,3))+1);
        else q(q_data(i,1)) = q(q_data(i,1)) + q_data(i,3);
        end
    else q(q_data(i,1)) = q(q_data(i,1)) + x(q_data(i,2))*q_data(i,3);
    end
end

for i = 1:size(q)
%     subs(q(i),x,ones(size(x)))
    if(abs(subs(q(i),x,zeros(size(x)))) < .000001)
        q(i) = x(k+1);
    end
end
% q

for i = 1:c_data_rows
    c(S_theta(i,1)) = c(S_theta(i,1)) + x(S_theta(i,2))*S_theta(i,3);
    if(i <= S_theta_rows) c(i) = c(i) - S_theta_RHS(i);
    end
end

G = [eye(h), -M, q];


% Solve Problem By Hand 
% Phase 1

%iteration 2
G_ = pivot(G,3,11);
G2 = pivot(G_,5,9);
temp = G2(3,:);
G2(3,:) = G2(5,:);
G2(5,:) = temp;
G2 = simplify(G2);

%iteration 3
G3 = simplify(pivot(G2,4,10));

%Phase 2
G1_2 = subs(G3, x(2), 0);

%iteration 2
G2_2 = simplify(pivot(G1_2,1,7));

%iteration 3
G3_2 = simplify(pivot(G1_2,6,12));

%iteration 4
G4_2 = simplify(pivot(G3_2,4,4));

%iteration 5
G5_2 = simplify(pivot(G4_2,3,3));

%iteration 6
G6_2 = simplify(pivot(G5_2,2,8));