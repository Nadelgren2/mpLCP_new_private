function regions = solve_mpLCP_5_5_17(h,k,M_data,q_data,S_theta,S_theta_RHS)

warning('off','MATLAB:rankDeficientMatrix');

t = clock;
max_time = 3600;
phase1_max_iter = 100;
phase2_max_iter = 100;
print_stuff = 0;
print_stuff2 = 0;
use_optimal_instead_of_feasible = 1;
simplify_g = 1;
plot_if_k_is_2 = 1;

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

G = [eye(h), -M, q]
%return
feas_pt = zeros(1,k);

% get a feasible point in S_theta

fun = @(x) 0;
x0 = zeros(1,k);
A = [];
b = [];
Aeq = [];
beq = [];
lb = zeros(1,k+2)-999999;
ub = zeros(1,k+2)+999999;
leq = matlabFunction(c,'Vars',x);
nonlcon = @(x) create_nonlcon(leq,[],x);
options = optimoptions('fmincon','Algorithm','active-set','Display', 'off','MaxIterations',4000,'MaxFunctionEvaluations',8000,'ConstraintTolerance',.000000001);
options_ = optimoptions('fmincon','Algorithm','active-set','MaxIterations',4000,'MaxFunctionEvaluations',8000,'ConstraintTolerance',.000000001);
options2 = optimoptions('fmincon','Algorithm','sqp','Display', 'off','ConstraintTolerance',.000000001);
options3 = optimoptions('fmincon','Algorithm','interior-point','Display', 'off','MaxIterations',4000,'MaxFunctionEvaluations',8000,'ConstraintTolerance',.000000001);
options2_ = optimoptions('fmincon','Algorithm','sqp','ConstraintTolerance',.000000001);
options3_ = optimoptions('fmincon','Algorithm','interior-point','MaxIterations',4000,'MaxFunctionEvaluations',8000,'ConstraintTolerance',.000000001);
[sol,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k),ub(1:k),nonlcon,options);

if(exitflag >= 1) feas_pt = [sol,1,0];
else printf('Feasible starting point not found');
    return
end

% Start phase 1

basis_stack = zeros(1,h+1); %the last element should indicate the row of discovered_bases, etc, which is associated with the basis
basis_stack(1,end) = 1;
phase1_iter = 1;
bases = basis_stack(1:h);
discovered_bases = basis_stack;
basis_dim = [1];
Z = zeros(1,h);
H = zeros(1,h^2);
E = zeros(1,h);
EQ = zeros(1,h);
F = zeros(1,h);
P = [1];
iota = [1];
kappa = [1];
% num_nonzero = [0];
RHS = sym(zeros(1,h));
pivots_made = [0,0,0];
pivot_from_prev = [0,0,0];
previous_basis = zeros(1,h+1);
previous_pivot = [0];
prev_fval = 1;
current_basis = [];
nlp_s_sol = [];
min_fval = 100000;
less_than_min_fval = 0;

% build Z,E,H for initial basis

%Z and create polynomial defining constraints of inv rgn
for i = 1:h
    indices(i) = i;
    [N,D] = numden(G(i,end));
    n = coeffs(N,'All');
    if(isempty(n))
        Z(1,i) = 1;
%         num_nonzero(1) = num_nonzero(1) + 1;
    else
        val = subs(D,x,feas_pt);
        if(size(n,1) > 1 || abs(n(1)) > .00001)
            if(val > 0) RHS(1,i) = -N;
            else RHS(1,i) = N;
            end
        end
    end
end

%E and H
for i = 1:h
    if(Z(1,i) ~= 1 && E(1,i) ~= 1 && size(coeffs(RHS(1,i),'All'),2) > 1)
        for j = 1:h
            if(Z(1,j) ~= 1 && E(1,j) ~= 1 && j ~= i)
                if(H(1,(i-1)*h + j) ~= 1)
                    fun = @(x) -x(k+2);
                    x0 = feas_pt(1:k+2);
                    temp = RHS(1,setdiff(indices,[i,j]));
                    leq = matlabFunction([c,temp(temp~=0),RHS(1,j)+x(k+2)],'Vars',x);
                    eq = matlabFunction(RHS(1,i),'Vars',x);
                    nonlcon = @(x) create_nonlcon(leq,eq,x);
                    [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                    if(exitflag == 0)
                        exitflag2 = 0;
                        [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                        if(exitflag2 ~= exitflag)
                            exitflag = exitflag2;
                        else
                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                            return
                        end
                    end
                    if(exitflag >= 1 && abs(fval) < .00001)
                        H(1,(i-1)*h + j) = 1;
                        for p = 1:h
                           if(H(1,(j-1)*h + p) == 1) 
                               H(1,(i-1)*h + p) = 1;
                           end
                        end
                    elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                        E(1,i) = 1;
                        break
                    end
                end
            end
        end
    end
end

%Phase 1 main step

phase2_feas_pt = [];

while(size(basis_stack,1) > 0 &&  phase1_iter < phase1_max_iter && etime(clock,t) < max_time)
   %initialize
   disp(sprintf('Phase 1 iteration: %d',phase1_iter));
   phase1_iter = phase1_iter + 1;
   current_basis = basis_stack(1,:);
   row = current_basis(1,end);
   basis_stack = basis_stack(2:end,:);
   
   %update G
   prev_G = G;
   if(print_stuff || print_stuff2)
       current_basis
      prev_G 
   end
   
   if(phase1_iter ~= 2)
       if(previous_pivot(row) == size(pivots_made,1))
           d_vs_e = pivot_from_prev(row,1);
           i = pivot_from_prev(row,2);
           j = pivot_from_prev(row,3);
           if(d_vs_e == 0)
               if(previous_basis(i) == 0)
                   G = pivot(G,i,i+h);
               else
                   G = pivot(G,i,i);
               end
               pivots_made = [pivots_made;0,i,0];
           else
               if(previous_basis(i) == 0)
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   end
               else
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   end
               end
               pivots_made = [pivots_made;1,i,j];
           end
       else
           while(previous_pivot(row) < size(pivots_made,1))
               undo_pivot = pivots_made(end,:);
               pivots_made = pivots_made(1:end-1,:);
               d_vs_e = undo_pivot(1);
               i = undo_pivot(2);
               j = undo_pivot(3);
               if(print_stuff)
                   d_vs_e
                   undo_pivot
                   i
                   j
               end
               if(d_vs_e == 0)
                   if(previous_basis(i) == 0)
                       G = pivot(G,i,i+h);
                       previous_basis(i) = 1;
                   else
                       G = pivot(G,i,i);
                       previous_basis(i) = 0;
                   end
               else
                   if(previous_basis(i) == 0)
                       if(previous_basis(j) == 0)
                           G = pivot(G,i,j+h);
                           G = pivot(G,j,i+h);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 1;
                           previous_basis(j) = 1;
                       else
                           G = pivot(G,i,j);
                           G = pivot(G,j,i+h);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 1;
                           previous_basis(j) = 0;
                       end
                   else
                       if(previous_basis(j) == 0)
                           G = pivot(G,i,j+h);
                           G = pivot(G,j,i);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 0;
                           previous_basis(j) = 1;
                       else
                           G = pivot(G,i,j);
                           G = pivot(G,j,i);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 0;
                           previous_basis(j) = 0;
                       end
                   end
               end
           end
           d_vs_e = pivot_from_prev(row,1);
           i = pivot_from_prev(row,2);
           j = pivot_from_prev(row,3);
           if(d_vs_e == 0)
               if(previous_basis(i) == 0)
                   G = pivot(G,i,i+h);
               else
                   G = pivot(G,i,i);
               end
               pivots_made = [pivots_made;0,i,0];
           else
               if(previous_basis(i) == 0)
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   end
               else
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   end
               end
               pivots_made = [pivots_made;1,i,j];
           end
       end
   end
   
   if(simplify_g)
       G = simplify(G);
   end
   
   if(print_stuff)
       G
   end
   
   
   %set up nlp_s.
   
%    c_new = c;
%    for i = 1:h
%        [N,D] = numden(G(i,end));
%        val = subs(D,x,feas_pt);
%        n = coeffs(N,'All');
%        if(~isempty(n))
%            if(size(n,1) > 1 || abs(n(1)) > .00001)
%                if(val > 0) c_new = [c_new,-N*D];
%                else c_new = [c_new,N*D];
%                end
%            end
%        else
%            Z(row,i) = 1;
%        end
%    end

%    basis_dim(row)
%    fval
   if(basis_dim(row) == 1)
       fun = @(x) x(k+1);
       x0 = feas_pt(1:k+2);
       if(length(phase2_feas_pt) > 0)
           x0 = phase2_feas_pt(1:k+2);
       end
       temp = RHS(row,:);
       leq = matlabFunction([c,temp(temp~=0)],'Vars',x);
       eq = [];
       nonlcon = @(x) create_nonlcon(leq,eq,x);
       if(print_stuff)
           x0
           temp(temp~=0)
           [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options_)
           [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2_)
           [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options3_)
       end
       [sol,fval,exitflag] = fmincon_4_nlp_s(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options,options2,options3);
       nlp_s_sol = sol;
       phase2_feas_pt = sol;
       %feas_pt = sol;
       if( abs(fval - min_fval) < .00001 )
           less_than_min_fval = 1;
       else
           less_than_min_fval = 0;
       end
       if(exitflag ~= 0)
           min_fval = min([min_fval,fval]);
           prev_fval = fval;
       end
   else
       fval = prev_fval;
       exitflag = 1;
   end
%    fval
   
   if(exitflag == 0)
       [sol,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
       nlp_s_sol = sol;
       if(exitflag2 ~= exitflag)
           exitflag = exitflag2;
           min_fval = min([min_fval,fval]);
           prev_fval = fval;
       end
   end
   
   disp(sprintf('Value of rho at nlp_s solution: %.10f',fval));
   if(exitflag < 1) 
%        if(print_stuff)
%            exitflag
%        end
       disp(sprintf('Fmincon failed. Exitflag: %d Line %d',exitflag, MFileLineNr()));
       G = prev_G;
   elseif(fval <= 0.000001)
       if(print_stuff)
            G(:,end)
       end
       if(basis_dim(row) == 1)
           %set up nlp_d (for phase 2)
           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
           x0 = feas_pt(1:k+2);
           temp = subs(RHS(row,:),x(k+1),0);
           leq = matlabFunction([c,temp(temp~=0)+x(k+2)],'Vars',x);
           eq = [];
           nonlcon = @(x) create_nonlcon(leq,eq,x);
           [sol,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
           if(exitflag >= 1)
               bases = current_basis(1:h);
               regions = [regions; temp];
%                phase2_feas_pt = sol;
               break
           else
               for i = 1:h
                  if(Z(row,i) ~= 1 && E(row,i) ~= 1 && F(row,i) ~= 1)
                     discard = [];
                     for j = 1:h
                        if(Z(row,j) == 1 || H(row,(i-1)*h + j) == 1 || j == i)
                            discard = [discard,j];
                        end
                     end
                     fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                     x0 = feas_pt(1:k+2);
                     temp = RHS(row,setdiff(indices,discard));
                     leq = matlabFunction([c,temp+x(k+2)],'Vars',x);
                     eq = matlabFunction(RHS(row,i),'Vars',x);
                     nonlcon = @(x) create_nonlcon(leq,eq,x);
                     [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                     if(exitflag >= 1) 
                         F(row,i) = 1;
                         for j = 1:h
                            if(H(row,(i-1)*h + j) == 1)
                               F(row,j) = 1; 
                            end
                         end
                     end
                  end
               end
               %find adjacent regions across each facet
               for l=1:h
                  if(F(row,l) == 1)
                      if(current_basis(l) == 0) % l is a w
                         if(G(l,l+h) ~= 0)
                            new_basis = current_basis;
                            new_basis(l) = 1;
                            if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                               new_basis(h+1) = size(Z,1)+1;
                               if(print_stuff)
                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                    new_basis
                               end
                               basis_stack = [new_basis; basis_stack];
                               discovered_bases = [new_basis; discovered_bases]; %more is needed
                               G_ = pivot(G,l,l+h);
                               temp_Z = zeros(1,h);
                               temp_E = zeros(1,h);
                               temp_H = zeros(1,h^2);
                               temp_RHS = sym(zeros(1,h));
                               %build Z and RHS
                               for m = 1:h
                                    [N,d] = numden(G_(m,end));
                                    n = coeffs(N,'All');
                                    if(isempty(n))
                                        temp_Z(1,m) = 1;
                                    else
                                        val = subs(d,x,feas_pt);
                                        if(size(n,1) > 1 || abs(n(1)) > .00001)
                                            if(val > 0) 
                                                temp_RHS(1,m) = -N;
                                            else
                                                temp_RHS(1,m) = N;
                                            end
                                        end
                                    end
                               end
                               %build E and H
                               for m = 1:h
                                     if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                        for n = 1:h
                                            if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                if(temp_H(1,(m-1)*h + n) ~= 1)
                                                    fun = @(x) -x(k+2);
                                                    x0 = feas_pt(1:k+2);
                                                    temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                    leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                    eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                    nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                    [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                    if(exitflag == 0)
                                                        [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                        if(exitflag2 ~= exitflag)
                                                            exitflag = exitflag2;
%                                                         else
%                                                             disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
%                                                             return
                                                        end
                                                    end
                                                    if(exitflag >=0 && abs(fval) < .00001)
                                                        temp_H(1,(m-1)*h + n) = 1;
                                                        for p = 1:h
                                                           if(temp_H(1,(n-1)*h + p) == 1) 
                                                               temp_H(1,(m-1)*h + p) = 1;
                                                           end
                                                        end
                                                    elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                        temp_E(1,m) = 1;
                                                        break
                                                    end
                                                end
                                             end
                                        end
                                     end
                               end %end of build E and H
                               Z = [Z; temp_Z];
                               E = [E; temp_E];
                               H = [H; temp_H];
                               EQ = [EQ; zeros(1,h)];
                               F = [F; zeros(1,h)];
                               RHS = [RHS; temp_RHS];
                               kappa = [kappa; l];
                               bases = [bases;new_basis(1:h)];
                               pivot_from_prev = [pivot_from_prev; 0,l,0];
                               previous_pivot = [previous_pivot; size(pivots_made,1)];
                               % build P and iota
                               if(basis_dim(row) == 1)
                                    P = [P;row];
                                    iota = [iota;l];
                               else
                                    P = [P;P(row)];
                                    iota = [iota;iota(row)];
                               end
                               % set up nlp_d
                               fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                               x0 = feas_pt(1:k+2);
                               leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                               nonlcon = @(x) create_nonlcon(leq,[],x);
                               if(print_stuff)
                                                nonlcon 
                                             end
                               [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                               if(exitflag >= 1)
                                     basis_dim = [basis_dim;1];
                               elseif(exitflag == -2)
                                     basis_dim = [basis_dim;0];
                               else
                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                    return  
                               end
                            end
                         else
                            for j = 1:h
                               if(j ~= l)
                                   if(current_basis(j) == 0) % l and j are both w's
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [1,1];
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l+h);
                                             G_ = pivot(G_,l,j+h);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) temp_RHS(1,m) = -N;
                                                        else temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+2);
                                                                x0 = feas_pt(1:k+2);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                             leq = matlabFunction([c,temp],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                if(print_stuff)
                                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                    new_basis
                                                end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 kappa = [kappa; j];
                                                 bases = [bases;new_basis(1:h)];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                                 x0 = feas_pt(1:k+2);
                                                 leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                                 if(print_stuff)
                                                nonlcon 
                                             end
                                                 if(exitflag >= 1)
                                                     basis_dim = [basis_dim;1];
                                                 elseif(exitflag == -2)
                                                     basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                                 else
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end
                                   else % l is a w, j is a z
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [1,0];
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l+h);
                                             G_ = pivot(G_,l,j);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) temp_RHS(1,m) = -N;
                                                        else temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+2);
                                                                x0 = feas_pt(1:k+2);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                             leq = matlabFunction([c,temp],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                 if(print_stuff)
                                                        disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                        new_basis
                                                 end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 kappa = [kappa; j];
                                                 bases = [bases;new_basis(1:h)];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                                 x0 = feas_pt(1:k+2);
                                                 leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 if(print_stuff)
                                                nonlcon 
                                             end
                                                 [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                                 if(exitflag >= 1)
                                                     basis_dim = [basis_dim;1];
                                                 elseif(exitflag == -2)
                                                     basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                                 else
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end
                                   end
                               end
                            end
                         end
                      else % l is a z
                          if(G(l,l) ~= 0)
                            new_basis = current_basis;
                            new_basis(l) = 0;
                            if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                               new_basis(h+1) = size(Z,1)+1;
                               if(print_stuff)
                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                    new_basis
                               end
                               basis_stack = [new_basis; basis_stack];
                               discovered_bases = [new_basis; discovered_bases]; %more is needed
                               G_ = pivot(G,l,l);
                               temp_Z = zeros(1,h);
                               temp_E = zeros(1,h);
                               temp_H = zeros(1,h^2);
                               temp_RHS = sym(zeros(1,h));
                               %build Z and RHS
                               for m = 1:h
                                    [N,d] = numden(G_(m,end));
                                    n = coeffs(N,'All');
                                    if(isempty(n))
                                        temp_Z(1,m) = 1;
                                    else
                                        val = subs(d,x,feas_pt);
                                        if(size(n,1) > 1 || abs(n(1)) > .00001)
                                            if(val > 0) 
                                                temp_RHS(1,m) = -N;
                                            else
                                                temp_RHS(1,m) = N;
                                            end
                                        end
                                    end
                               end
                               %build E and H
                               for m = 1:h
                                     if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                        for n = 1:h
                                            if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                if(temp_H(1,(m-1)*h + n) ~= 1)
                                                    fun = @(x) -x(k+2);
                                                    x0 = feas_pt(1:k+2);
                                                    temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                    leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                    eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                    nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                    [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                    if(exitflag == 0)
                                                        exitflag2 = 0;
                                                        [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                        if(exitflag2 ~= exitflag)
                                                            exitflag = exitflag2;
                                                        else
                                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                            return
                                                        end
                                                    end
                                                    if(exitflag >= 1 && abs(fval) < .00001)
                                                        temp_H(1,(m-1)*h + n) = 1;
                                                        for p = 1:h
                                                           if(temp_H(1,(n-1)*h + p) == 1) 
                                                               temp_H(1,(m-1)*h + p) = 1;
                                                           end
                                                        end
                                                    elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                        temp_E(1,m) = 1;
                                                        break
                                                    end
                                                end
                                             end
                                        end
                                     end
                               end %end of build E and H
                               Z = [Z; temp_Z];
                               E = [E; temp_E];
                               H = [H; temp_H];
                               EQ = [EQ; zeros(1,h)];
                               F = [F; zeros(1,h)];
                               RHS = [RHS; temp_RHS];
                               kappa = [kappa; l];
                               bases = [bases;new_basis(1:h)];
                               pivot_from_prev = [pivot_from_prev; 0,l,0];
                               previous_pivot = [previous_pivot; size(pivots_made,1)];
                               % build P and iota
                               if(basis_dim(row) == 1)
                                    P = [P;row];
                                    iota = [iota;l];
                               else
                                    P = [P;P(row)];
                                    iota = [iota;iota(row)];
                               end
                               % set up nlp_d
                               fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                               x0 = feas_pt(1:k+2);
                               leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                               nonlcon = @(x) create_nonlcon(leq,[],x);
                               if(print_stuff)
                                                nonlcon 
                                             end
                               [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                               if(exitflag >= 1)
                                     basis_dim = [basis_dim;1];
                               elseif(exitflag == -2)
                                     basis_dim = [basis_dim;0];
                               else
                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                    return  
                               end
                            end
                         else
                            for j = 1:h
                               if(j ~= l)
                                   if(current_basis(j) == 0) % l is a z, and j is a w
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [0,1];
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l);
                                             G_ = pivot(G_,l,j+h);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) temp_RHS(1,m) = -N;
                                                        else temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+2);
                                                                x0 = feas_pt(1:k+2);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                             leq = matlabFunction([c,temp],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                 if(print_stuff)
                                                        disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                        new_basis
                                                 end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 kappa = [kappa; j];
                                                 bases = [bases;new_basis(1:h)];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                                 x0 = feas_pt(1:k+2);
                                                 leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 if(print_stuff)
                                                nonlcon 
                                             end
                                                 [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                                 if(exitflag >= 1)
                                                     basis_dim = [basis_dim;1];
                                                 elseif(exitflag == -2)
                                                     basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                                 else
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end
                                   else % l is a z, j is a z
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [0,0];
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l);
                                             G_ = pivot(G_,l,j);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) temp_RHS(1,m) = -N;
                                                        else temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+2);
                                                                x0 = feas_pt(1:k+2);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                             leq = matlabFunction([c,temp],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                 if(print_stuff)
                                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                    new_basis
                                                 end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 kappa = [kappa; j];
                                                 bases = [bases;new_basis(1:h)];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                                 x0 = feas_pt(1:k+2);
                                                 leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 if(print_stuff)
                                                nonlcon 
                                             end
                                                 [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                                 if(exitflag >= 1)
                                                     basis_dim = [basis_dim;1];
                                                 elseif(exitflag == -2)
                                                     basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                                 else
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end
                                   end
                               end
                            end
                         end
                      end
                  end
               end %end find adjacent
           end
       end
   else
       if(print_stuff)
            disp(sprintf('Reached Line %d',MFileLineNr()));
            G(:,end)
       end
       %basis_dim_row = basis_dim(row)
       if(basis_dim(row) == 1)
           %find phase 1 facets
           for i = 1:h
               if(abs(subs(RHS(row,i),x,sol)) < .00001)
                   EQ(row,i) = 1;
               end
           end
           for i = 1:h
              if(EQ(row,i) == 1 && Z(row,i) ~= 1 && E(row,i) ~= 1 && F(row,i) ~= 1)
                 discard = [];
                 for j = 1:h
                    if(Z(row,j) == 1 || H(row,(i-1)*h + j) == 1 || j == i)
                        discard = [discard,j];
                    end
                 end
                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                 x0 = feas_pt(1:k+2);
                 temp = RHS(row,setdiff(indices,discard));
                 leq = matlabFunction([c,temp+x(k+2)],'Vars',x);
                 eq = matlabFunction(RHS(row,i),'Vars',x);
                 nonlcon = @(x) create_nonlcon(leq,eq,x);
                 [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                 if(exitflag >= 1) 
                     F(row,i) = 1;
                     for j = 1:h
                        if(H(row,(i-1)*h + j) == 1)
                           F(row,j) = 1; 
                        end
                     end
                 end
              end
           end
           if(print_stuff)
               disp(sprintf('Line %d',MFileLineNr()));
                F
           end
           %sort F elements based on gradient
           sorted_facets = [];
           if(sum(F(row,:)) > 1)
               for i = 1:h
                  if(F(row,i) == 1)
                     if(print_stuff)
                         disp(sprintf('Line %d',MFileLineNr()));
                         subs(gradient(RHS(row,i),x),x,nlp_s_sol)
                     end
                     g = double(subs(gradient(RHS(row,i),x),x,nlp_s_sol));
                     g = g/norm(g);
                     sorted_facets = [sorted_facets;i,g(k+1)];
                  end
               end
               sorted_facets = sortrows(sorted_facets,2);
           else
               for i = 1:h
                  if(F(row,i) == 1)
                      sorted_facets = [i,0];
                  end
               end
           end
           if(print_stuff)
               disp(sprintf('Line %d',MFileLineNr()));
              sorted_facets 
              
           end
           %find adjacent regions across each facet
           for i = size(sorted_facets,1):-1:1
              l = sorted_facets(i,1);
              if(print_stuff)
                  l
              end
              if(current_basis(l) == 0) % l is a w
                  if(print_stuff)
                     G(l,l+h) 
                  end
                 if(G(l,l+h) ~= 0)
                    new_basis = current_basis;
                    new_basis(l) = 1;
                    if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                       new_basis(h+1) = size(Z,1)+1;
                       if(print_stuff)
                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                            new_basis
                       end
                       basis_stack = [new_basis; basis_stack];
                       discovered_bases = [new_basis; discovered_bases]; %more is needed
                       G_ = pivot(G,l,l+h);
                       temp_Z = zeros(1,h);
                       temp_E = zeros(1,h);
                       temp_H = zeros(1,h^2);
                       temp_RHS = sym(zeros(1,h));
                       %build Z and RHS
                       for m = 1:h
                            [N,d] = numden(G_(m,end));
                            n = coeffs(N,'All');
                            if(isempty(n))
                                temp_Z(1,m) = 1;
                            else
                                val = subs(d,x,feas_pt);
                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                    if(val > 0) 
                                        temp_RHS(1,m) = -N;
                                    else
                                        temp_RHS(1,m) = N;
                                    end
                                end
                            end
                       end
                       %build E and H
                       for m = 1:h
                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                for n = 1:h
                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                            fun = @(x) -x(k+2);
                                            x0 = feas_pt(1:k+2);
                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
%                                             if(exitflag == 0)
%                                                 exitflag2 = 0;
%                                                 [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
%                                                 if(exitflag2 ~= exitflag)
%                                                     exitflag = exitflag2;
%                                                 elseif(exitflag2 == 0)
%                                                     exitflag3 = 0;
%                                                     options4 = optimoptions('fmincon','Algorithm','interior-point','Display', 'off','MaxIterations',5000,'MaxFunctionEvaluations',10000);
%                                                     [~,fval,exitflag3] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options4);
%                                                     if(exitflag3 ~= exitflag)
%                                                         exitflag = exitflag3;
%                                                     else
%                                                         disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag3,MFileLineNr()));
%                                                         return
%                                                     end
%                                                 end
%                                             end
                                            if(exitflag >= 0 && abs(fval) < .00001)
                                                temp_H(1,(m-1)*h + n) = 1;
                                                for p = 1:h
                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                       temp_H(1,(m-1)*h + p) = 1;
                                                   end
                                                end
                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                temp_E(1,m) = 1;
                                                break
                                            end
                                        end
                                     end
                                end
                             end
                       end %end of build E and H
                       Z = [Z; temp_Z];
                       E = [E; temp_E];
                       H = [H; temp_H];
                       EQ = [EQ; zeros(1,h)];
                       F = [F; zeros(1,h)];
                       RHS = [RHS; temp_RHS];
                       kappa = [kappa; l];
                       bases = [bases;new_basis(1:h)];
                       pivot_from_prev = [pivot_from_prev; 0,l,0];
                       previous_pivot = [previous_pivot; size(pivots_made,1)];
                       % build P and iota
                       if(basis_dim(row) == 1)
                            P = [P;row];
                            iota = [iota;l];
                       else
                            P = [P;P(row)];
                            iota = [iota;iota(row)];
                       end
                       % set up nlp_d
                       fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                       x0 = feas_pt(1:k+2);
                       leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                       nonlcon = @(x) create_nonlcon(leq,[],x);
                       if(print_stuff)
                                                nonlcon 
                                             end
                       [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                       if(exitflag >= 1)
                             basis_dim = [basis_dim;1];
                       elseif(exitflag == -2)
                             basis_dim = [basis_dim;0];
                       else
                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                            return  
                       end
                    end
                 else
                    for j = 1:h
                       if(j ~= l)
                           if(print_stuff)
                               disp(sprintf('Line %d',MFileLineNr()));
                               j
                           end
                           if(current_basis(j) == 0) % l and j are both w's
                              new_basis = current_basis;
                              new_basis([l,j]) = [1,1];
                              if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                  if(print_stuff)
                                      disp(sprintf('Line %d',MFileLineNr()));
                                      subs(G(j,l+h),x,feas_pt)
                                      subs(G(l,j+h),x,feas_pt)
                                  end
                                  if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                     G_ = pivot(G,j,l+h);
                                     G_ = pivot(G_,l,j+h);
                                     G_([l,j],:) = G_([j,l],:);
                                     temp_Z = zeros(1,h);
                                     temp_E = zeros(1,h);
                                     temp_H = zeros(1,h^2);
                                     temp_RHS = sym(zeros(1,h));
                                     %build Z and RHS
                                     for m = 1:h
                                        [N,d] = numden(G_(m,end));
                                        n = coeffs(N,'All');
                                        if(isempty(n))
                                            temp_Z(1,m) = 1;
                                        else
                                            val = subs(d,x,feas_pt);
                                            if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                if(val > 0) temp_RHS(1,m) = -N;
                                                else temp_RHS(1,m) = N;
                                                end
                                            end
                                        end
                                     end
                                     %build E and H
                                     for m = 1:h
                                         if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                            for n = 1:h
                                                if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                    if(temp_H(1,(m-1)*h + n) ~= 1)
                                                        fun = @(x) -x(k+2);
                                                        x0 = feas_pt(1:k+2);
                                                        temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                        leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                        eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                        nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                        [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                        if(exitflag == 0)
                                                            exitflag2 = 0;
                                                            [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                            if(exitflag2 ~= exitflag)
                                                                exitflag = exitflag2;
                                                            else
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                        end
                                                        if(exitflag >= 1 && abs(fval) < .00001)
                                                            temp_H(1,(m-1)*h + n) = 1;
                                                            for p = 1:h
                                                               if(temp_H(1,(n-1)*h + p) == 1) 
                                                                   temp_H(1,(m-1)*h + p) = 1;
                                                               end
                                                            end
                                                        elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                            temp_E(1,m) = 1;
                                                            break
                                                        end
                                                    end
                                                 end
                                            end
                                         end
                                     end %end of build E and H
                                     %set up nlp_a
                                     discard1 = [];
                                     discard2 = [];
                                     for m = 1:h
                                        if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                            discard1 = [discard1,m];
                                        end
                                        if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                            discard2 = [discard2,m];
                                        end
                                     end
                                     fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                     x0 = feas_pt(1:k+2);
                                     temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                     leq = matlabFunction([c,temp],'Vars',x);
                                     eq = matlabFunction(RHS(row,l),'Vars',x);
                                     nonlcon = @(x) create_nonlcon(leq,eq,x);
                                     [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options,options2,options3);
                                     if(exitflag == 0)
                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                        return
                                     elseif(exitflag >= 1)
                                         new_basis(h+1) = size(Z,1)+1;
                                         if(print_stuff)
                                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                            new_basis
                                         end
                                         basis_stack = [new_basis; basis_stack];
                                         discovered_bases = [new_basis; discovered_bases];
                                         Z = [Z; temp_Z];
                                         E = [E; temp_E];
                                         H = [H; temp_H];
                                         EQ = [EQ; zeros(1,h)];
                                         F = [F; zeros(1,h)];
                                         RHS = [RHS; temp_RHS];
                                         kappa = [kappa; j];
                                         bases = [bases;new_basis(1:h)];
                                         pivot_from_prev = [pivot_from_prev; 1,l,j];
                                         previous_pivot = [previous_pivot; size(pivots_made,1)];
                                         % build P and iota
                                         if(basis_dim(row) == 1)
                                            P = [P;row];
                                            iota = [iota;l];
                                         else
                                            P = [P;P(row)];
                                            iota = [iota;iota(row)];
                                         end
                                         % set up nlp_d
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,[],x);
                                         if(print_stuff)
                                                nonlcon 
                                             end
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag >= 1)
                                             basis_dim = [basis_dim;1];
                                         elseif(exitflag == -2)
                                             basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                         else
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return  
                                         end
                                     end                                     
                                  end
                              end
                           else % l is a w, j is a z
                              new_basis = current_basis;
                              new_basis([l,j]) = [1,0];
                              if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                  if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                     G_ = pivot(G,j,l+h);
                                     G_ = pivot(G_,l,j);
                                     G_([l,j],:) = G_([j,l],:);
                                     temp_Z = zeros(1,h);
                                     temp_E = zeros(1,h);
                                     temp_H = zeros(1,h^2);
                                     temp_RHS = sym(zeros(1,h));
                                     %build Z and RHS
                                     for m = 1:h
                                        [N,d] = numden(G_(m,end));
                                        n = coeffs(N,'All');
                                        if(isempty(n))
                                            temp_Z(1,m) = 1;
                                        else
                                            val = subs(d,x,feas_pt);
                                            if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                if(val > 0) temp_RHS(1,m) = -N;
                                                else temp_RHS(1,m) = N;
                                                end
                                            end
                                        end
                                     end
                                     %build E and H
                                     for m = 1:h
                                         if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                            for n = 1:h
                                                if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                    if(temp_H(1,(m-1)*h + n) ~= 1)
                                                        fun = @(x) -x(k+2);
                                                        x0 = feas_pt(1:k+2);
                                                        temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                        leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                        eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                        nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                        [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                        if(exitflag == 0)
                                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                            return
                                                        end
                                                        if(exitflag >= 1 && abs(fval) < .00001)
                                                            temp_H(1,(m-1)*h + n) = 1;
                                                            for p = 1:h
                                                               if(temp_H(1,(n-1)*h + p) == 1) 
                                                                   temp_H(1,(m-1)*h + p) = 1;
                                                               end
                                                            end
                                                        elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                            temp_E(1,m) = 1;
                                                            break
                                                        end
                                                    end
                                                 end
                                            end
                                         end
                                     end %end of build E and H
                                     %set up nlp_a
                                     discard1 = [];
                                     discard2 = [];
                                     for m = 1:h
                                        if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                            discard1 = [discard1,m];
                                        end
                                        if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                            discard2 = [discard2,m];
                                        end
                                     end
                                     fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                     x0 = feas_pt(1:k+2);
                                     temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                     leq = matlabFunction([c,temp],'Vars',x);
                                     eq = matlabFunction(RHS(row,l),'Vars',x);
                                     nonlcon = @(x) create_nonlcon(leq,eq,x);
                                     [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                     if(exitflag == 0)
                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                        return
                                     elseif(exitflag >= 1)
                                         new_basis(h+1) = size(Z,1)+1;
                                         if(print_stuff)
                                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                            new_basis
                                         end
                                         basis_stack = [new_basis; basis_stack];
                                         discovered_bases = [new_basis; discovered_bases];
                                         Z = [Z; temp_Z];
                                         E = [E; temp_E];
                                         H = [H; temp_H];
                                         EQ = [EQ; zeros(1,h)];
                                         F = [F; zeros(1,h)];
                                         RHS = [RHS; temp_RHS];
                                         kappa = [kappa; j];
                                         bases = [bases;new_basis(1:h)];
                                         pivot_from_prev = [pivot_from_prev; 1,l,j];
                                         previous_pivot = [previous_pivot; size(pivots_made,1)];
                                         % build P and iota
                                         if(basis_dim(row) == 1)
                                            P = [P;row];
                                            iota = [iota;l];
                                         else
                                            P = [P;P(row)];
                                            iota = [iota;iota(row)];
                                         end
                                         % set up nlp_d
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,[],x);
                                         if(print_stuff)
                                                nonlcon 
                                             end
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag >= 1)
                                             basis_dim = [basis_dim;1];
                                         elseif(exitflag == -2)
                                             basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                         else
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return  
                                         end
                                     end                                     
                                  end
                              end
                           end
                       end
                    end
                 end
              else % l is a z
                  if(G(l,l) ~= 0)
                    new_basis = current_basis;
                    new_basis(l) = 0;
                    if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                       new_basis(h+1) = size(Z,1)+1;
                       if(print_stuff)
                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                            new_basis
                       end
                       basis_stack = [new_basis; basis_stack];
                       discovered_bases = [new_basis; discovered_bases]; %more is needed
                       G_ = pivot(G,l,l);
                       temp_Z = zeros(1,h);
                       temp_E = zeros(1,h);
                       temp_H = zeros(1,h^2);
                       temp_RHS = sym(zeros(1,h));
                       %build Z and RHS
                       for m = 1:h
                            [N,d] = numden(G_(m,end));
                            n = coeffs(N,'All');
                            if(isempty(n))
                                temp_Z(1,m) = 1;
                            else
                                val = subs(d,x,feas_pt);
                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                    if(val > 0) 
                                        temp_RHS(1,m) = -N;
                                    else
                                        temp_RHS(1,m) = N;
                                    end
                                end
                            end
                       end
                       %build E and H
                       for m = 1:h
                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                for n = 1:h
                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                            fun = @(x) -x(k+2);
                                            x0 = feas_pt(1:k+2);
                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                            if(exitflag == 0)
                                                exitflag2 = 0;
                                                [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                if(exitflag2 ~= exitflag)
                                                    exitflag = exitflag2;
                                                else
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return
                                                end
                                            end
                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                temp_H(1,(m-1)*h + n) = 1;
                                                for p = 1:h
                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                       temp_H(1,(m-1)*h + p) = 1;
                                                   end
                                                end
                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                temp_E(1,m) = 1;
                                                break
                                            end
                                        end
                                     end
                                end
                             end
                       end %end of build E and H
                       Z = [Z; temp_Z];
                       E = [E; temp_E];
                       H = [H; temp_H];
                       EQ = [EQ; zeros(1,h)];
                       F = [F; zeros(1,h)];
                       RHS = [RHS; temp_RHS];
                       kappa = [kappa; l];
                       bases = [bases;new_basis(1:h)];
                       pivot_from_prev = [pivot_from_prev; 0,l,0];
                       previous_pivot = [previous_pivot; size(pivots_made,1)];
                       % build P and iota
                       if(basis_dim(row) == 1)
                            P = [P;row];
                            iota = [iota;l];
                       else
                            P = [P;P(row)];
                            iota = [iota;iota(row)];
                       end
                       % set up nlp_d
                       fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                       x0 = feas_pt(1:k+2);
                       leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                       nonlcon = @(x) create_nonlcon(leq,[],x);
                       if(print_stuff)
                                                nonlcon 
                                             end
                       [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                       if(exitflag >= 1)
                             basis_dim = [basis_dim;1];
                       elseif(exitflag == -2)
                             basis_dim = [basis_dim;0];
                       else
                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                            return  
                       end
                    end
                 else
                    for j = 1:h
                       if(j ~= l)
                           if(current_basis(j) == 0) % l is a z, and j is a w
                              new_basis = current_basis;
                              new_basis([l,j]) = [0,1];
                              if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                  if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                     G_ = pivot(G,j,l);
                                     G_ = pivot(G_,l,j+h);
                                     G_([l,j],:) = G_([j,l],:);
                                     temp_Z = zeros(1,h);
                                     temp_E = zeros(1,h);
                                     temp_H = zeros(1,h^2);
                                     temp_RHS = sym(zeros(1,h));
                                     %build Z and RHS
                                     for m = 1:h
                                        [N,d] = numden(G_(m,end));
                                        n = coeffs(N,'All');
                                        if(isempty(n))
                                            temp_Z(1,m) = 1;
                                        else
                                            val = subs(d,x,feas_pt);
                                            if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                if(val > 0) temp_RHS(1,m) = -N;
                                                else temp_RHS(1,m) = N;
                                                end
                                            end
                                        end
                                     end
                                     %build E and H
                                     for m = 1:h
                                         if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                            for n = 1:h
                                                if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                    if(temp_H(1,(m-1)*h + n) ~= 1)
                                                        fun = @(x) -x(k+2);
                                                        x0 = feas_pt(1:k+2);
                                                        temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                        leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                        eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                        nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                        [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                        if(exitflag == 0)
                                                            exitflag2 = 0;
                                                            [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                            if(exitflag2 ~= exitflag)
                                                                exitflag = exitflag2;
                                                            else
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                        end
                                                        if(exitflag >= 1 && abs(fval) < .00001)
                                                            temp_H(1,(m-1)*h + n) = 1;
                                                            for p = 1:h
                                                               if(temp_H(1,(n-1)*h + p) == 1) 
                                                                   temp_H(1,(m-1)*h + p) = 1;
                                                               end
                                                            end
                                                        elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                            temp_E(1,m) = 1;
                                                            break
                                                        end
                                                    end
                                                 end
                                            end
                                         end
                                     end %end of build E and H
                                     %set up nlp_a
                                     discard1 = [];
                                     discard2 = [];
                                     for m = 1:h
                                        if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                            discard1 = [discard1,m];
                                        end
                                        if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                            discard2 = [discard2,m];
                                        end
                                     end
                                     fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                     x0 = feas_pt(1:k+2);
                                     temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                     leq = matlabFunction([c,temp],'Vars',x);
                                     eq = matlabFunction(RHS(row,l),'Vars',x);
                                     nonlcon = @(x) create_nonlcon(leq,eq,x);
                                     [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                     if(exitflag == 0)
                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                        return
                                     elseif(exitflag >= 1)
                                         new_basis(h+1) = size(Z,1)+1;
                                         if(print_stuff)
                                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                            new_basis
                                         end
                                         basis_stack = [new_basis; basis_stack];
                                         discovered_bases = [new_basis; discovered_bases];
                                         Z = [Z; temp_Z];
                                         E = [E; temp_E];
                                         H = [H; temp_H];
                                         EQ = [EQ; zeros(1,h)];
                                         F = [F; zeros(1,h)];
                                         RHS = [RHS; temp_RHS];
                                         kappa = [kappa; j];
                                         bases = [bases;new_basis(1:h)];
                                         pivot_from_prev = [pivot_from_prev; 1,l,j];
                                         previous_pivot = [previous_pivot; size(pivots_made,1)];
                                         % build P and iota
                                         if(basis_dim(row) == 1)
                                            P = [P;row];
                                            iota = [iota;l];
                                         else
                                            P = [P;P(row)];
                                            iota = [iota;iota(row)];
                                         end
                                         % set up nlp_d
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,[],x);
                                         if(print_stuff)
                                                nonlcon 
                                             end
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag >= 1)
                                             basis_dim = [basis_dim;1];
                                         elseif(exitflag == -2)
                                             basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                         else
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return  
                                         end
                                     end                                     
                                  end
                              end
                           else % l is a z, j is a z
                              new_basis = current_basis;
                              new_basis([l,j]) = [0,0];
                              if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                  if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                     G_ = pivot(G,j,l);
                                     G_ = pivot(G_,l,j);
                                     G_([l,j],:) = G_([j,l],:);
                                     temp_Z = zeros(1,h);
                                     temp_E = zeros(1,h);
                                     temp_H = zeros(1,h^2);
                                     temp_RHS = sym(zeros(1,h));
                                     %build Z and RHS
                                     for m = 1:h
                                        [N,d] = numden(G_(m,end));
                                        n = coeffs(N,'All');
                                        if(isempty(n))
                                            temp_Z(1,m) = 1;
                                        else
                                            val = subs(d,x,feas_pt);
                                            if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                if(val > 0) temp_RHS(1,m) = -N;
                                                else temp_RHS(1,m) = N;
                                                end
                                            end
                                        end
                                     end
                                     %build E and H
                                     for m = 1:h
                                         if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                            for n = 1:h
                                                if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                    if(temp_H(1,(m-1)*h + n) ~= 1)
                                                        fun = @(x) -x(k+2);
                                                        x0 = feas_pt(1:k+2);
                                                        temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                        leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                        eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                        nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                        [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                        if(exitflag == 0)
                                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                            return
                                                        end
                                                        if(exitflag >= 1 && abs(fval) < .00001)
                                                            temp_H(1,(m-1)*h + n) = 1;
                                                            for p = 1:h
                                                               if(temp_H(1,(n-1)*h + p) == 1) 
                                                                   temp_H(1,(m-1)*h + p) = 1;
                                                               end
                                                            end
                                                        elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                            temp_E(1,m) = 1;
                                                            break
                                                        end
                                                    end
                                                 end
                                            end
                                         end
                                     end %end of build E and H
                                     %set up nlp_a
                                     discard1 = [];
                                     discard2 = [];
                                     for m = 1:h
                                        if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                            discard1 = [discard1,m];
                                        end
                                        if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                            discard2 = [discard2,m];
                                        end
                                     end
                                     fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                     x0 = feas_pt(1:k+2);
                                     temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                     leq = matlabFunction([c,temp],'Vars',x);
                                     eq = matlabFunction(RHS(row,l),'Vars',x);
                                     nonlcon = @(x) create_nonlcon(leq,eq,x);
                                     [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                     if(exitflag == 0)
                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                        return
                                     elseif(exitflag >= 1)
                                         new_basis(h+1) = size(Z,1)+1;
                                         if(print_stuff)
                                            disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                            new_basis
                                         end
                                         basis_stack = [new_basis; basis_stack];
                                         discovered_bases = [new_basis; discovered_bases];
                                         Z = [Z; temp_Z];
                                         E = [E; temp_E];
                                         H = [H; temp_H];
                                         EQ = [EQ; zeros(1,h)];
                                         F = [F; zeros(1,h)];
                                         RHS = [RHS; temp_RHS];
                                         kappa = [kappa; j];
                                         bases = [bases;new_basis(1:h)];
                                         pivot_from_prev = [pivot_from_prev; 1,l,j];
                                         previous_pivot = [previous_pivot; size(pivots_made,1)];
                                         % build P and iota
                                         if(basis_dim(row) == 1)
                                            P = [P;row];
                                            iota = [iota;l];
                                         else
                                            P = [P;P(row)];
                                            iota = [iota;iota(row)];
                                         end
                                         % set up nlp_d
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,[],x);
                                         if(print_stuff)
                                                nonlcon 
                                         end
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag >= 1)
                                             basis_dim = [basis_dim;1];
                                         elseif(exitflag == -2)
                                             basis_dim = [basis_dim;0]; %#ok<*AGROW>
                                         else
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return  
                                         end
                                     end                                     
                                  end
                              end
                           end
                       end
                    end
                 end
              end
           end %end find adjacent
           %determine number new
           number_new = size(Z,1) - row;
           if(number_new == 0 && less_than_min_fval == 1 && min_fval > 0)
               disp(sprintf('The Problem is Infeasible, exitting!!')); %#ok<*DSPS>
               return
           end
           %eliminate overlapping adjacent
           for i = 1:h
              if(Z(row,i) == 1)
                 if(current_basis(i) == 0)
                     if(G(i,i+h) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 1;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 else
                     if(G(i,i) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 0;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 end
              end
           end
       else % current inv rgn is low dim
           % get GCDs
           if(print_stuff)
               disp(sprintf('Current inv region is low dim'));
           end
           GCD = sym(zeros(1,h));
           D = zeros(1,h);
           for i = 1:h
               if(i ~= kappa(row))
                   GCD(i) = gcd(RHS(row,i),RHS(row,kappa(row)));
               end
           end
           %GCD
           %D
           for i = 1:h
              n = coeffs(GCD(i),'All'); 
              %size_ = size(n,1)
              if(~isempty(n) && (size(n,1) > 1 || abs(n(1)) > .00001))
                 %set up nlp_g
                 discard = [];
                 for m = 1:h
                    if(Z(row,m) == 1 || m == i)
                        discard = [discard,m];
                    end
                 end
                 fun = @(x) -x(k+2);
                 x0 = feas_pt(1:k+2);
                 temp = RHS(row,setdiff(indices,discard));
                 leq = matlabFunction([c,temp,GCD(i)^2+x(k+2)],'Vars',x);
                 eq = matlabFunction(RHS(row,i),'Vars',x);
                 nonlcon = @(x) create_nonlcon(leq,eq,x);
                 [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                 if(print_stuff)
                     disp(sprintf('Line %d',MFileLineNr()));
                    i
                    sol
                    fval
                 end
                 if(exitflag == 0)
                     [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                     if(exitflag == 0)
                         disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                         return
                     end
                 end
                 if(exitflag >= 1 && abs(fval) < .00001)
                    D(i) = 1;
                 end
              end
           end  
           if(print_stuff)
             disp(sprintf('Line %d',MFileLineNr()));
             D
           end
           %find adjacent regions across each facet in D
           for l = 1:h
              if(D(l) == 1) 
                  if(current_basis(l) == 0) % l is a w
                     if(G(l,l+h) ~= 0)
                        new_basis = current_basis;
                        new_basis(l) = 1;
                        if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                           new_basis(h+1) = size(Z,1)+1;
                           if(print_stuff)
                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                new_basis
                           end
                           basis_stack = [new_basis; basis_stack];
                           discovered_bases = [new_basis; discovered_bases]; 
                           G_ = pivot(G,l,l+h);
                           temp_Z = zeros(1,h);
                           temp_E = zeros(1,h);
                           temp_H = zeros(1,h^2);
                           temp_RHS = sym(zeros(1,h));
                           %build Z and RHS
                           for m = 1:h
                                [N,d] = numden(G_(m,end));
                                n = coeffs(N,'All');
                                if(isempty(n))
                                    temp_Z(1,m) = 1;
                                else
                                    val = subs(d,x,feas_pt);
                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                        if(val > 0) 
                                            temp_RHS(1,m) = -N;
                                        else
                                            temp_RHS(1,m) = N;
                                        end
                                    end
                                end
                           end
                           %build E and H
                           for m = 1:h
                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                    for n = 1:h
                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                fun = @(x) -x(k+2);
                                                x0 = feas_pt(1:k+2);
                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                if(exitflag == 0)
                                                    exitflag2 = 0;
                                                    [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                    if(exitflag2 ~= exitflag)
                                                        exitflag = exitflag2;
                                                    else
                                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                        return
                                                    end
                                                end
                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                    temp_H(1,(m-1)*h + n) = 1;
                                                    for p = 1:h
                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                           temp_H(1,(m-1)*h + p) = 1;
                                                       end
                                                    end
                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                    temp_E(1,m) = 1;
                                                    break
                                                end
                                            end
                                         end
                                    end
                                 end
                           end %end of build E and H
                           Z = [Z; temp_Z];
                           E = [E; temp_E];
                           H = [H; temp_H];
                           EQ = [EQ; zeros(1,h)];
                           F = [F; zeros(1,h)];
                           RHS = [RHS; temp_RHS];
                           kappa = [kappa; l];
                           bases = [bases;new_basis(1:h)];
                           pivot_from_prev = [pivot_from_prev; 0,l,0];
                           previous_pivot = [previous_pivot; size(pivots_made,1)];
                           % build P and iota
                           if(basis_dim(row) == 1)
                                P = [P;row];
                                iota = [iota;l];
                           else
                                P = [P;P(row)];
                                iota = [iota;iota(row)];
                           end
                           % set up nlp_d
                           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                           x0 = feas_pt(1:k+2);
                           leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                           nonlcon = @(x) create_nonlcon(leq,[],x);
                           if(print_stuff)
                            nonlcon 
                           end
                           [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                           if(exitflag >= 1)
                                 basis_dim = [basis_dim;1];
                           elseif(exitflag == -2)
                                 basis_dim = [basis_dim;0];
                           else
                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                return  
                           end
                        end
                     else
                        for j = 1:h
                           if(j ~= l)
                               if(current_basis(j) == 0) % l and j are both w's
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [1,1];
                                  if(print_stuff)
                                      disp(sprintf('Checking basis:'))
                                      new_basis
                                      disp(sprintf('Line %d',MFileLineNr()));
                                      is_in(new_basis(1:h),discovered_bases(:,1:h))
                                      if(current_basis(end) == 4)
                                      discovered_bases(13,:)
                                      end
                                      subs(G(j,l+h),x,feas_pt)
                                      subs(G(l,j+h),x,feas_pt)
                                  end
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j+h);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+2);
                                                            x0 = feas_pt(1:k+2);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                         leq = matlabFunction([c,temp],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         if(print_stuff)
                                             [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options)
                                         else
                                             [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options,options2,options3);
                                         end
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             kappa = [kappa; j];
                                             bases = [bases;new_basis(1:h)];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             if(print_stuff)
                                                nonlcon 
                                             end
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag >= 1)
                                                 basis_dim = [basis_dim;1];
                                             elseif(exitflag == -2)
                                                 basis_dim = [basis_dim;0];
                                             else
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               else % l is a w, j is a z
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [1,0];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+2);
                                                            x0 = feas_pt(1:k+2);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                            if(exitflag == 0)
                                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                         leq = matlabFunction([c,temp],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             kappa = [kappa; j];
                                             bases = [bases;new_basis(1:h)];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag >= 1)
                                                 basis_dim = [basis_dim;1];
                                             elseif(exitflag == -2)
                                                 basis_dim = [basis_dim;0];
                                             else
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end 
                               end
                           end
                        end
                     end
                  else % l is a z
                      if(G(l,l) ~= 0)
                        new_basis = current_basis;
                        new_basis(l) = 0;
                        if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                           new_basis(h+1) = size(Z,1)+1;
                           if(print_stuff)
                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                new_basis
                           end
                           basis_stack = [new_basis; basis_stack];
                           discovered_bases = [new_basis; discovered_bases] ;
                           G_ = pivot(G,l,l);
                           temp_Z = zeros(1,h);
                           temp_E = zeros(1,h);
                           temp_H = zeros(1,h^2);
                           temp_RHS = sym(zeros(1,h));
                           %build Z and RHS
                           for m = 1:h
                                [N,d] = numden(G_(m,end));
                                n = coeffs(N,'All');
                                if(isempty(n))
                                    temp_Z(1,m) = 1;
                                else
                                    val = subs(d,x,feas_pt);
                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                        if(val > 0) 
                                            temp_RHS(1,m) = -N;
                                        else
                                            temp_RHS(1,m) = N;
                                        end
                                    end
                                end
                           end
                           %build E and H
                           for m = 1:h
                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                    for n = 1:h
                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                fun = @(x) -x(k+2);
                                                x0 = feas_pt(1:k+2);
                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                if(exitflag == 0)
                                                    [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options2);
                                                    if(exitflag == 0)
                                                        disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                        return
                                                    end
                                                end
                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                    temp_H(1,(m-1)*h + n) = 1;
                                                    for p = 1:h
                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                           temp_H(1,(m-1)*h + p) = 1;
                                                       end
                                                    end
                                                elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                    temp_E(1,m) = 1;
                                                    break
                                                end
                                            end
                                         end
                                    end
                                 end
                           end %end of build E and H
                           Z = [Z; temp_Z];
                           E = [E; temp_E];
                           H = [H; temp_H];
                           EQ = [EQ; zeros(1,h)];
                           F = [F; zeros(1,h)];
                           RHS = [RHS; temp_RHS];
                           kappa = [kappa; l];
                           bases = [bases;new_basis(1:h)];
                           pivot_from_prev = [pivot_from_prev; 0,l,0];
                           previous_pivot = [previous_pivot; size(pivots_made,1)];
                           % build P and iota
                           if(basis_dim(row) == 1)
                                P = [P;row];
                                iota = [iota;l];
                           else
                                P = [P;P(row)];
                                iota = [iota;iota(row)];
                           end
                           % set up nlp_d
                           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                           x0 = feas_pt(1:k+2);
                           leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                           nonlcon = @(x) create_nonlcon(leq,[],x);
                           [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                           if(exitflag >= 1)
                                 basis_dim = [basis_dim;1];
                           elseif(exitflag == -2)
                                 basis_dim = [basis_dim;0];
                           else
                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                return  
                           end
                        end
                     else
                        for j = 1:h
                           if(j ~= l)
                               if(current_basis(j) == 0) % l is a z and j is a w
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [0,1];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l);
                                         G_ = pivot(G_,l,j+h);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+2);
                                                            x0 = feas_pt(1:k+2);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                         leq = matlabFunction([c,temp],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             kappa = [kappa; j];
                                             bases = [bases;new_basis(1:h)];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             leq = matlabFunction([c,temp_RHS(temp_RHS ~= 0) + x(k+2)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             if(print_stuff)
                                                 temp_RHS
                                                nonlcon 
                                                [sol,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options_)
                                             end
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag >= 1)
                                                 basis_dim = [basis_dim;1];
                                             elseif(exitflag == -2)
                                                 basis_dim = [basis_dim;0];
                                             else
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               else % l is a z, j is a z
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [0,0];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l);
                                         G_ = pivot(G_,l,j);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+2);
                                                            x0 = feas_pt(1:k+2);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+2)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+2),nonlcon,options);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                         x0 = feas_pt(1:k+2);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+2);
                                         leq = matlabFunction([c,temp],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             kappa = [kappa; j];
                                             bases = [bases;new_basis(1:h)];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
                                             x0 = feas_pt(1:k+2);
                                             leq = matlabFunction([c,temp_RHS + x(k+2)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [~,~,exitflag] = fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k+1),.0001],ub(1:k+2),nonlcon,options);
                                             if(exitflag >= 1)
                                                 basis_dim = [basis_dim;1];
                                             elseif(exitflag == -2)
                                                 basis_dim = [basis_dim;0];
                                             else
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               end
                           end
                        end
                     end
                  end
              end
           end %end find adjacent
           %determine number new
%            number_new = size(Z,1) - row;
%            if(number_new == 0)
%                disp(sprintf('The Problem is Infeasible, exitting!!'));
%                return
%            end
           %eliminate overlapping adjacent
           for i = 1:h
              if(Z(row,i) == 1)
                 if(current_basis(i) == 0)
                     if(G(i,i+h) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 1;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 else
                     if(G(i,i) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 0;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 end
              end
           end
       end
   end
   %colors = [];
   %G(:,end);
   %regions_ = G(:,end);
   previous_basis = current_basis;
   %figure;
   %hold on;
   %colors = plot_RHS(regions_,c,h,x(1:2),colors);
   %return
end

 if(min_fval > 0)
       disp(sprintf('The Problem is Infeasible, exitting!!'));
       return
 end

% disp(sprintf('Phase 1 iteration: %d',phase1_iter));

colors = [];

if(k == 2 && plot_if_k_is_2)
    figure;
    hold on;
    colors = plot_RHS(regions,c,h,x,colors);
end

if(print_stuff)
    c
end

regions = simplify(subs(G(:,end),[x(k+1),x(k+2)],[0,0]));
phase2_feas_pt = [phase2_feas_pt(1:k),phase2_feas_pt(k+2)];
bases = current_basis(1:h);

if(phase1_iter ~= 2 && 0)
    if(previous_pivot(row) == size(pivots_made,1))
        d_vs_e = pivot_from_prev(row,1);
        i = pivot_from_prev(row,2);
        j = pivot_from_prev(row,3);
        if(print_stuff)
            d_vs_e
            i
            j
            G
        end
        if(d_vs_e == 0)
           if(previous_basis(i) == 0)
               G = pivot(G,i,i+h);
           else
               G = pivot(G,i,i);
           end
        else
           if(previous_basis(i) == 0)
               if(previous_basis(j) == 0)
                   G = pivot(G,i,j+h);
                   G = pivot(G,j,i+h);
                   G([i,j],:) = G([j,i],:);
               else
                   G = pivot(G,i,j);
                   G = pivot(G,j,i+h);
                   G([i,j],:) = G([j,i],:);
               end
           else
               if(previous_basis(j) == 0)
                   G = pivot(G,i,j+h);
                   G = pivot(G,j,i);
                   G([i,j],:) = G([j,i],:);
               else
                   G = pivot(G,i,j);
                   G = pivot(G,j,i);
                   G([i,j],:) = G([j,i],:);
               end
           end
        end
    else
        disp(sprintf('doh. Line %d',MFileLineNr()));
        return
    end
end

y = sym('y',[1 k+1]);
G = subs(G,x,[y(1:k),0,y(k+1)]);
x = sym('x',[1 k+1]);
G = subs(G,y,x);

current_basis(end) = 1;
basis_stack = current_basis;
discovered_bases = basis_stack;
basis_dim = [1]; %#ok<*NBRAK>
Z = zeros(1,h);
H = zeros(1,h^2);
E = zeros(1,h);
F = zeros(1,h);
P = [1];
iota = [1];
kappa = [1];
RHS = sym(zeros(1,h));
pivots_made = [0,0,0];
pivot_from_prev = [0,0,0];
previous_basis = zeros(1,h+1);
previous_pivot = [0];
run_phase_2 = 1;
phase2_iter = 1;
feas_pt = [feas_pt(1:k),feas_pt(k+2)];
lb = [lb(1:k),lb(k+2)];
ub = [ub(1:k),ub(k+2)];

% build Z and E and H for start of phase 2
for i = 1:h
    indices(i) = i;
    [N,D] = numden(G(i,end));
    n = coeffs(N,'All');
    if(isempty(n))
        Z(1,i) = 1;
    else
        val = subs(D,x,feas_pt);
        if(size(n,1) > 1 || abs(n(1)) > .00001)
            if(val > 0) 
                RHS(1,i) = -N;
            else
                RHS(1,i) = N;
            end
        end
    end
end

%try to get an interior point
if(min(abs(phase2_feas_pt(1:k))) < .01)
    fun = @(x) -x(k+1);
    x0 = phase2_feas_pt;
    leq = matlabFunction([c+x(k+1),RHS(1,:)+x(k+1)],'Vars',x);
    eq = [];
    nonlcon = @(x) create_nonlcon(leq,eq,x);
    if(print_stuff)
        [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3)
    else
        [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
    end
    if(exitflag >= 1)
        phase2_feas_pt = sol;
    end
end

nlp_d_pts = zeros(h,k+1);
for i = 1:h
    nlp_d_pts(i,:) = phase2_feas_pt;
end

%E and H
for i = 1:h
    if(print_stuff)
        i
    end
%     Z(1,i)
%     E(1,i)
%     coeffs(RHS(1,i),'All')
%     size(coeffs(RHS(1,i),'All'),1)
    if(Z(1,i) ~= 1 && E(1,i) ~= 1 && size(coeffs(RHS(1,i),'All'),2) > 1)
        for j = 1:h
            if(print_stuff)
                 j
            end
            if(Z(1,j) ~= 1 && E(1,j) ~= 1 && j ~= i)
                if(H(1,(i-1)*h + j) ~= 1)
                    fun = @(x) -x(k+1);
                    x0 = phase2_feas_pt;
                    temp = RHS(1,setdiff(indices,[i,j]));
                    leq = matlabFunction([c+x(k+1),temp(temp~=0),RHS(1,j)+x(k+1)],'Vars',x);
                    eq = matlabFunction(RHS(1,i),'Vars',x);
                    nonlcon = @(x) create_nonlcon(leq,eq,x);
                    if(print_stuff)
                        [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3)
                    else
                        [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                    end
%                     if(exitflag == 0)
%                         exitflag2 = 0;
%                         [~,fval,exitflag2] = fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options2);
%                         if(exitflag2 ~= exitflag)
%                             exitflag = exitflag2;
%                         else
%                             disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
%                             return
%                         end
%                     end
                    if(exitflag >= 1 && abs(fval) < .00001)
                        H(1,(i-1)*h + j) = 1;
                        for p = 1:h
                           if(H(1,(j-1)*h + p) == 1) 
                               H(1,(i-1)*h + p) = 1;
                           end
                        end
                    elseif(exitflag == -2 || (exitflag >= 1 && -fval < 0))
                        E(1,i) = 1;
                        break
                    end
                end
            end
        end
    end
end
if(print_stuff)
    disp(Z); %#ok<*UNRCH>
    disp(E);
    disp(H);
end

sol = zeros(1,k+1);
pivots_made = [0,0,0];
pivot_from_prev = [0,0,0];
feas_pts = x0;

% begin phase 2
while(size(basis_stack,1) > 0 &&  phase2_iter < phase2_max_iter && etime(clock,t) < max_time)
   %initialize
   disp(sprintf('Phase 2 iteration: %d. Elapsed time: %.2f',phase2_iter,etime(clock,t)));
   phase2_iter = phase2_iter + 1;
   current_basis = basis_stack(1,:);
   row = current_basis(1,end);
   basis_stack = basis_stack(2:end,:);
   x0 = feas_pts(row,:);
   
   if(print_stuff)
       current_basis
       row
   end
   
   %update G
%    prev_G = G;
   
   if(phase2_iter ~= 2)
%        row
%        previous_pivot
%        pivots_made
%        size(pivots_made,1)
       if(previous_pivot(row) == size(pivots_made,1))
           d_vs_e = pivot_from_prev(row,1);
           i = pivot_from_prev(row,2);
           j = pivot_from_prev(row,3);
           if(d_vs_e == 0)
               if(previous_basis(i) == 0)
                   G = pivot(G,i,i+h);
               else
                   G = pivot(G,i,i);
               end
               pivots_made = [pivots_made;0,i,0];
           else
               if(previous_basis(i) == 0)
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   end
               else
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   end
               end
               pivots_made = [pivots_made;1,i,j];
           end
       else
           while(previous_pivot(row) < size(pivots_made,1))
               undo_pivot = pivots_made(end,:);
               pivots_made = pivots_made(1:end-1,:);
               d_vs_e = undo_pivot(1);
               i = undo_pivot(2);
               j = undo_pivot(3);
               if(print_stuff)
                   d_vs_e
                   undo_pivot
                   i
                   j
               end
               if(d_vs_e == 0)
                   if(previous_basis(i) == 0)
                       G = pivot(G,i,i+h);
                       previous_basis(i) = 1;
                   else
                       G = pivot(G,i,i);
                       previous_basis(i) = 0;
                   end
               else
                   if(previous_basis(i) == 0)
                       if(previous_basis(j) == 0)
                           G = pivot(G,i,j+h);
                           G = pivot(G,j,i+h);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 1;
                           previous_basis(j) = 1;
                       else
                           G = pivot(G,i,j);
                           G = pivot(G,j,i+h);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 1;
                           previous_basis(j) = 0;
                       end
                   else
                       if(previous_basis(j) == 0)
                           G = pivot(G,i,j+h);
                           G = pivot(G,j,i);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 0;
                           previous_basis(j) = 1;
                       else
                           G = pivot(G,i,j);
                           G = pivot(G,j,i);
                           G([i,j],:) = G([j,i],:);
                           previous_basis(i) = 0;
                           previous_basis(j) = 0;
                       end
                   end
               end
           end
           d_vs_e = pivot_from_prev(row,1);
           i = pivot_from_prev(row,2);
           j = pivot_from_prev(row,3);
           if(d_vs_e == 0)
               if(previous_basis(i) == 0)
                   G = pivot(G,i,i+h);
               else
                   G = pivot(G,i,i);
               end
               pivots_made = [pivots_made;0,i,0];
           else
               if(previous_basis(i) == 0)
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i+h);
                       G([i,j],:) = G([j,i],:);
                   end
               else
                   if(previous_basis(j) == 0)
                       G = pivot(G,i,j+h);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   else
                       G = pivot(G,i,j);
                       G = pivot(G,j,i);
                       G([i,j],:) = G([j,i],:);
                   end
               end
               pivots_made = [pivots_made;1,i,j];
           end
       end
   end
   if(print_stuff)
     simplify(G)
     simplify(G(:,end))
   end
   
   if(simplify_g)
       G = simplify(G);
   end
   
   %try to get an interior point
%     fun = @(x) -x(k+1);
%     leq = matlabFunction([c+x(k+1),RHS(row,:)+x(k+1)],'Vars',x);
%     eq = [];
%     nonlcon = @(x) create_nonlcon(leq,eq,x);
%     if(print_stuff)
%         [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3)
%     else
%         [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
%     end
%     if(exitflag >= 1)
%         if(basis_dim(row) ~= 1)
%             basis_dim(row) = 1;
%         end
%         phase2_feas_pt = sol;
%     else
%         phase2_feas_pt = x0;
%     end
  
   if( run_phase_2 )
%        current_basis
%        basis_dim
       if(basis_dim(row) == 1)
           %find phase 2 facets
           if(print_stuff)
               disp('Full dim.');
           end
           for i = 1:h
               if(print_stuff)
                    i
               end
              if(Z(row,i) ~= 1 && E(row,i) ~= 1 && F(row,i) ~= 1)
                 discard = [];
                 for j = 1:h
                    if(Z(row,j) == 1 || H(row,(i-1)*h + j) == 1 || j == i)
                        discard = [discard,j];
                    end
                 end
                 %nlp_f
                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                  x0 = feas_pt(1:k+1);
%                  x0
                 temp = RHS(row,setdiff(indices,discard));
                 leq = matlabFunction([c+x(k+1),temp(temp~=0)+x(k+1)],'Vars',x);
                 eq = matlabFunction(RHS(row,i),'Vars',x);
                 nonlcon = @(x) create_nonlcon(leq,eq,x);
                 exitflag = 0;
                 if(print_stuff)
                     [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3)
                 else
                    [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                 end
                 if(exitflag >= 1) 
                     nlp_d_pts(i,:) = sol;
                     F(row,i) = 1;
                     for j = 1:h
                        if(H(row,(i-1)*h + j) == 1)
                           F(row,j) = 1; 
                        end
                     end
                 end
              end
           end
           if(print_stuff)
               disp(Z);
               disp(E);
               disp(H);
               disp(F);
           end
           %find adjacent regions across each facet
           for l = 1:h
%                l
%                F(row,l)
              if(F(row,l) == 1)
%                   current_basis(l)
                  if(current_basis(l) == 0) % l is a w
                     if(print_stuff)
                           disp('l is a w.');
                     end
                     if(G(l,l+h) ~= 0)
                         if(print_stuff)
                           disp('diagonal pivot.');
                         end
                        new_basis = current_basis;
                        new_basis(l) = 1;
                        if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                           new_basis(h+1) = size(Z,1)+1;
                           if(print_stuff)
                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                new_basis
                           end
                           basis_stack = [new_basis; basis_stack];
                           discovered_bases = [new_basis; discovered_bases];
                           G_ = pivot(G,l,l+h);
                           temp_Z = zeros(1,h);
                           temp_E = zeros(1,h);
                           temp_H = zeros(1,h^2);
                           temp_RHS = sym(zeros(1,h));
%                            new_pt = x0;
                           %build Z and RHS
                           for m = 1:h
                                [N,d] = numden(G_(m,end));
                                if(print_stuff)
                                    G_(m,end)
                                    N
                                    d
                                end
                                n = coeffs(N,'All');
                                if(isempty(n))
                                    temp_Z(1,m) = 1;
                                else
                                    val = subs(d,x,feas_pt);
                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                        if(val > 0) 
                                            temp_RHS(1,m) = -N;
                                        else
                                            temp_RHS(1,m) = N;
                                        end
                                    end
                                end
                           end
                           %build E and H
                           for m = 1:h
                               if(print_stuff) 
                                   m 
                               end
                               one_was_feas = 0;
                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                    for n = 1:h
                                        if(print_stuff)
                                            n
                                        end
                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                fun = @(x) -x(k+1);
%                                                 x0
%                                                 x0 = feas_pt(1:k+1);
%                                                 x0 = [.4,.05,0]
                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                if(print_stuff)
                                                    [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3)
                                                else
                                                    [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                end
                                                if(exitflag >= 1 && -fval >= 0)
                                                    one_was_feas = 1;
                                                end
                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                    temp_H(1,(m-1)*h + n) = 1;
                                                    for p = 1:h
                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                           temp_H(1,(m-1)*h + p) = 1;
                                                       end
                                                    end
                                                elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                    temp_E(1,m) = 1;
                                                    break
                                                end
                                            end
                                         end
                                    end
                                 end
                           end %end of build E and H
                           Z = [Z; temp_Z];
                           E = [E; temp_E];
                           H = [H; temp_H];
                           F = [F; zeros(1,h)];
                           RHS = [RHS; temp_RHS];
                           if(print_stuff)
                              temp_Z
                              temp_E
                              temp_H
                              temp_RHS
                           end
                           if(k == 2 && plot_if_k_is_2)
                             colors = plot_RHS(temp_RHS,c,h,x,colors);
                           end
                           kappa = [kappa; l];
                           pivot_from_prev = [pivot_from_prev; 0,l,0];
                           previous_pivot = [previous_pivot; size(pivots_made,1)];
                           % build P and iota
                           if(basis_dim(row) == 1)
                                P = [P;row];
                                iota = [iota;l];
                           else
                                P = [P;P(row)];
                                iota = [iota;iota(row)];
                           end
                           % set up nlp_d
                           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                            x0 = feas_pt(1:k+1);
                           leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                           nonlcon = @(x) create_nonlcon(leq,[],x);
                           exitflag = 0;
                           if(print_stuff)
                               temp_RHS
                               leq
                                nlp_d_pts(l,:)
                                [sol,fval,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.00001],ub(1:k+1),nonlcon,options,options2,options3)
                                %disp(sprintf('Sol: %.8f\t %.8f\t %.8f',sol(1),sol(2),sol(3)));
                           else
                                [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.00001],ub(1:k+1),nonlcon,options,options2,options3);
                           end
                           if(exitflag >= 1)
                                 feas_pts = [feas_pts;sol];
                                 basis_dim = [basis_dim;1];
                                 bases = [bases;new_basis(1:h)];
                                 regions = [regions;simplify(G_(:,end))];
                           elseif(exitflag == -2)
                                 feas_pts = [feas_pts;x0];
                                 basis_dim = [basis_dim;0];
                           else
                                 feas_pts = [feas_pts;x0];
                                 basis_dim = [basis_dim;0];
%                                 disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
%                                 return  
                           end
                        end
                     else
                         if(print_stuff)
                           disp('check excahnge pivots.');
                         end
                        for j = 1:h
%                             j
                           if(j ~= l)
                               if(current_basis(j) == 0) % l and j are both w's
                                 if(print_stuff)
                                       disp('l and j are both w_s.');
                                 end
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [1,1];
%                                   is_in(new_basis(1:h),discovered_bases(:,1:h))
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
%                                       G
%                                       subs(G(j,l+h),x,feas_pt)
%                                       subs(G(l,j+h),x,feas_pt)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j+h);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif( one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0]; 
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               else % l is a w, j is a z
                                 if(print_stuff)
                                       disp('l is a w, j is a z.');
                                 end
                                 new_basis = current_basis;
                                  new_basis([l,j]) = [1,0];
%                                   is_in(new_basis(1:h),discovered_bases(:,1:h))
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
%                                       G
%                                       subs(G(j,l+h),x,feas_pt)
%                                       subs(G(l,j+h),x,feas_pt)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0]; 
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               end
                           end
                        end
                     end
                  else % l is a z
                      if(current_basis(l) == 1) % l is a z
                          if(print_stuff)
                           disp('l is a z.');
                          end
                         if(print_stuff)
                             G(l,l)
                         end
                         if(G(l,l) ~= 0)
                            new_basis = current_basis;
                            new_basis(l) = 0;
                            if(print_stuff)
                               is_in(new_basis(1:h),discovered_bases(:,1:h))
                               discovered_bases
                            end
                            if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                               new_basis(h+1) = size(Z,1)+1;
                               if(print_stuff)
                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                    new_basis
                               end
                               basis_stack = [new_basis; basis_stack];
                               discovered_bases = [new_basis; discovered_bases];
                               G_ = pivot(G,l,l);
                               temp_Z = zeros(1,h);
                               temp_E = zeros(1,h);
                               temp_H = zeros(1,h^2);
                               temp_RHS = sym(zeros(1,h));
%                                new_pt = x0;
                               %build Z and RHS
                               for m = 1:h
                                    [N,d] = numden(G_(m,end));
                                    n = coeffs(N,'All');
                                    if(isempty(n))
                                        temp_Z(1,m) = 1;
                                    else
                                        val = subs(d,x,feas_pt);
                                        if(size(n,1) > 1 || abs(n(1)) > .00001)
                                            if(val > 0) 
                                                temp_RHS(1,m) = -N;
                                            else
                                                temp_RHS(1,m) = N;
                                            end
                                        end
                                    end
                               end
                               %build E and H
                               for m = 1:h
                                   one_was_feas = 0;
                                   if(print_stuff)
                                       m
                                   end
                                     if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                        for n = 1:h
                                            if(print_stuff)
                                                n
                                            end
                                            if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                if(temp_H(1,(m-1)*h + n) ~= 1)
                                                    fun = @(x) -x(k+1);
%                                                     x0 = feas_pt(1:k+1);
                                                    temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                    leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                    eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                    nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                    if(print_stuff)
                                                        [sol,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),-1],ub(1:k+1),nonlcon,options,options2,options3)
                                                    else
                                                        [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                    end
%                                                     if(exitflag == 0)
%                                                         disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
%                                                         return
%                                                     end
                                                    if(exitflag >= 1 && -fval >= 0)
                                                        one_was_feas = 1;
                                                    end
                                                    if(exitflag >= 1 && abs(fval) < .00001)
                                                        temp_H(1,(m-1)*h + n) = 1;
                                                        for p = 1:h
                                                           if(temp_H(1,(n-1)*h + p) == 1) 
                                                               temp_H(1,(m-1)*h + p) = 1;
                                                           end
                                                        end
                                                    elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                        temp_E(1,m) = 1;
                                                        break
                                                    end
                                                end
                                             end
                                        end
                                     end
                               end %end of build E and H
                               Z = [Z; temp_Z];
                               E = [E; temp_E];
                               H = [H; temp_H];
                               EQ = [EQ; zeros(1,h)];
                               F = [F; zeros(1,h)];
                               RHS = [RHS; temp_RHS];
                               if(k == 2 && plot_if_k_is_2)
                                 colors = plot_RHS(temp_RHS,c,h,x,colors);
                               end
                               kappa = [kappa; l];
                               pivot_from_prev = [pivot_from_prev; 0,l,0];
                               previous_pivot = [previous_pivot; size(pivots_made,1)];
                               % build P and iota
                               if(basis_dim(row) == 1)
                                    P = [P;row];
                                    iota = [iota;l];
                               else
                                    P = [P;P(row)];
                                    iota = [iota;iota(row)];
                               end
                               % set up nlp_d
                               fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                x0 = feas_pt(1:k+1);
                               leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                               nonlcon = @(x) create_nonlcon(leq,[],x);
                               [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                               if(exitflag >= 1)
                                 feas_pts = [feas_pts;sol];
                                 basis_dim = [basis_dim;1];
                                 bases = [bases;new_basis(1:h)];
                                 regions = [regions;simplify(G_(:,end))];
                               elseif(exitflag == -2)
                                 feas_pts = [feas_pts;x0];
                                 basis_dim = [basis_dim;0];
                               else
                                feas_pts = [feas_pts;x0]; 
                                basis_dim = [basis_dim;0];
                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                return  
                               end
                            end
                         else
                            for j = 1:h
                                if(print_stuff)
                                    j
                                end
                               if(j ~= l)
                                   if(current_basis(j) == 0) % l is a z and j is a w
                                       if(print_stuff)
                                           disp('l is a z, j is a w.');
                                       end
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [0,1];
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if(print_stuff)
                                              subs(G(j,l),x,feas_pt) 
                                              subs(G(l,j+h),x,feas_pt)
                                          end
                                          if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l);
                                             G_ = pivot(G_,l,j+h);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
%                                              new_pt = x0;
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) 
                                                            temp_RHS(1,m) = -N;
                                                        else
                                                            temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 one_was_feas = 0;
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+1);
%                                                                 x0 = feas_pt(1:k+1);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && -fval >= 0)
                                                                    one_was_feas = 1;
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                             leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                 if(print_stuff)
                                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                    new_basis
                                                 end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 if(k == 2 && plot_if_k_is_2)
                                                    colors = plot_RHS(temp_RHS,c,h,x,colors);
                                                 end
                                                 kappa = [kappa; j];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                                  x0 = feas_pt(1:k+1);
                                                 leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                                 if(exitflag >= 1)
                                                     feas_pts = [feas_pts;sol];
                                                     basis_dim = [basis_dim;1];
                                                     bases = [bases;new_basis(1:h)];
                                                     regions = [regions;simplify(G_(:,end))];
                                                 elseif(exitflag == -2)
                                                     feas_pts = [feas_pts;x0];
                                                     basis_dim = [basis_dim;0];
                                                 else
                                                    feas_pts = [feas_pts;x0]; 
                                                    basis_dim = [basis_dim;0];
                                                    disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end
                                   else % l is a z, j is a z
                                      if(print_stuff)
                                           disp('l is a z, j is a z.');
                                       end
                                      new_basis = current_basis;
                                      new_basis([l,j]) = [0,0];
%                                       is_in(new_basis(1:h),discovered_bases(:,1:h))
                                      if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                          if(print_stuff)
                                              subs(G(j,l),x,feas_pt) 
                                              subs(G(l,j),x,feas_pt)
                                          end
                                          if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                             G_ = pivot(G,j,l);
                                             G_ = pivot(G_,l,j);
                                             G_([l,j],:) = G_([j,l],:);
                                             temp_Z = zeros(1,h);
                                             temp_E = zeros(1,h);
                                             temp_H = zeros(1,h^2);
                                             temp_RHS = sym(zeros(1,h));
%                                              new_pt = x0;
                                             %build Z and RHS
                                             for m = 1:h
                                                [N,d] = numden(G_(m,end));
                                                n = coeffs(N,'All');
                                                if(isempty(n))
                                                    temp_Z(1,m) = 1;
                                                else
                                                    val = subs(d,x,feas_pt);
                                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                        if(val > 0) 
                                                            temp_RHS(1,m) = -N;
                                                        else
                                                            temp_RHS(1,m) = N;
                                                        end
                                                    end
                                                end
                                             end
                                             %build E and H
                                             for m = 1:h
                                                 one_was_feas = 0;
                                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                    for n = 1:h
                                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                                fun = @(x) -x(k+1);
%                                                                 x0 = feas_pt(1:k+1);
                                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                                leq = matlabFunction([c+x(k+1),temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                                [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                                if(exitflag == 0)
                                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                    return
                                                                end
                                                                if(exitflag >= 1 && -fval >= 0)
                                                                    one_was_feas = 1;
                                                                end
                                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                                    temp_H(1,(m-1)*h + n) = 1;
                                                                    for p = 1:h
                                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                                           temp_H(1,(m-1)*h + p) = 1;
                                                                       end
                                                                    end
                                                                elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                    temp_E(1,m) = 1;
                                                                    break
                                                                end
                                                            end
                                                         end
                                                    end
                                                 end
                                             end %end of build E and H
                                             %set up nlp_a
                                             discard1 = [];
                                             discard2 = [];
                                             for m = 1:h
                                                if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                    discard1 = [discard1,m];
                                                end
                                                if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                    discard2 = [discard2,m];
                                                end
                                             end
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                             leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                             eq = matlabFunction(RHS(row,l),'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,eq,x);
                                             [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag == 0)
                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                return
                                             elseif(exitflag >= 1)
                                                 new_basis(h+1) = size(Z,1)+1;
                                                 if(print_stuff)
                                                    disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                    new_basis
                                                 end
                                                 basis_stack = [new_basis; basis_stack];
                                                 discovered_bases = [new_basis; discovered_bases];
                                                 Z = [Z; temp_Z];
                                                 E = [E; temp_E];
                                                 H = [H; temp_H];
                                                 EQ = [EQ; zeros(1,h)];
                                                 F = [F; zeros(1,h)];
                                                 RHS = [RHS; temp_RHS];
                                                 if(k == 2 && plot_if_k_is_2)
                                                    colors = plot_RHS(temp_RHS,c,h,x,colors);
                                                 end
                                                 kappa = [kappa; j];
                                                 pivot_from_prev = [pivot_from_prev; 1,l,j];
                                                 previous_pivot = [previous_pivot; size(pivots_made,1)];
                                                 % build P and iota
                                                 if(basis_dim(row) == 1)
                                                    P = [P;row];
                                                    iota = [iota;l];
                                                 else
                                                    P = [P;P(row)];
                                                    iota = [iota;iota(row)];
                                                 end
                                                 % set up nlp_d
                                                 fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                                  x0 = feas_pt(1:k+1);
                                                 leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                                 nonlcon = @(x) create_nonlcon(leq,[],x);
                                                 [sol,~,exitflag] = use_fmincon(fun,nlp_d_pts(l,:),A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                                 if(exitflag >= 1)
                                                     feas_pts = [feas_pts;sol];
                                                     basis_dim = [basis_dim;1];
                                                     bases = [bases;new_basis(1:h)];
                                                     regions = [regions;simplify(G_(:,end))];
                                                 elseif(exitflag == -2)
                                                     feas_pts = [feas_pts;x0];
                                                     basis_dim = [basis_dim;0];
                                                 else
                                                    feas_pts = [feas_pts;x0]; 
                                                    basis_dim = [basis_dim;0];
                                                    disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                    return  
                                                 end
                                             end                                     
                                          end
                                      end 
                                   end
                               end
                            end
                         end
                      end
                  end
              end
           end %end find adjacent
           %eliminate overlapping adjacent
           for i = 1:h
              if(Z(row,i) == 1)
                 if(current_basis(i) == 0)
                     if(G(i,i+h) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 1;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 else
                     if(G(i,i) ~= 0)
                        new_basis = current_basis;
                        new_basis(i) = 0;
                        discovered_bases = [new_basis; discovered_bases];
                     end
                 end
              end
           end
       else % current inv rgn is low dim
           if(print_stuff)
               disp('Low dim.');
           end
           % get GCDs
           GCD = sym(zeros(1,h));
           D = zeros(1,h);
           for i = 1:h
               if(i ~= kappa(row))
                   GCD(i) = gcd(RHS(row,i),RHS(row,kappa(row)));
               end
           end
           for i = 1:h
              n = coeffs(GCD(i),'All'); 
              if(size(n,1) > 1)
                 %set up nlp_g
                 discard = [];
                 for m = 1:h
                    if(Z(row,m) == 1 || m == i)
                        discard = [discard,m];
                    end
                 end
                 fun = @(x) -x(k+1);
%                  x0 = feas_pt(1:k+1);
                 temp = RHS(row,setdiff(indices,discard));
                 leq = matlabFunction([c,temp(temp~=0),GCD(i)^2+x(k+1)],'Vars',x);
                 eq = matlabFunction(RHS(row,i),'Vars',x);
                 nonlcon = @(x) create_nonlcon(leq,eq,x);
                 [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                 if(exitflag == 0)
                     disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                     return
                 end
                 if(exitflag >= 1 && abs(fval) < .00001)
                    D(i) = 1;
                 end
              end
           end  
           if(print_stuff)
               kappa
              disp(D);
              disp(GCD);
           end
           %find adjacent regions across each facet in D
           for l = 1:h
              if(D(l) == 1) 
                  if(current_basis(l) == 0) % l is a w
                     if(G(l,l+h) ~= 0)
                        new_basis = current_basis;
                        new_basis(l) = 1; 
                        if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                           new_basis(h+1) = size(Z,1)+1;
                           if(print_stuff)
                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                new_basis
                           end
                           basis_stack = [new_basis; basis_stack];
                           discovered_bases = [new_basis; discovered_bases]; 
                           G_ = pivot(G,l,l+h);
                           temp_Z = zeros(1,h);
                           temp_E = zeros(1,h);
                           temp_H = zeros(1,h^2);
                           temp_RHS = sym(zeros(1,h));
%                            new_pt = x0;
                           %build Z and RHS
                           for m = 1:h
                                [N,d] = numden(G_(m,end));
                                n = coeffs(N,'All');
                                if(isempty(n))
                                    temp_Z(1,m) = 1;
                                else
                                    val = subs(d,x,feas_pt);
                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                        if(val > 0) 
                                            temp_RHS(1,m) = -N;
                                        else
                                            temp_RHS(1,m) = N;
                                        end
                                    end
                                end
                           end
                           %build E and H
                           for m = 1:h
                               one_was_feas = 0;
                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                    for n = 1:h
                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                fun = @(x) -x(k+1);
%                                                 x0 = feas_pt(1:k+1);
                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
%                                                 if(exitflag == 0)
%                                                     disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
%                                                     return
%                                                 end
                                                if(exitflag >= 1 && -fval >= 0)
                                                    one_was_feas = 1;
                                                end
                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                    temp_H(1,(m-1)*h + n) = 1;
                                                    for p = 1:h
                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                           temp_H(1,(m-1)*h + p) = 1;
                                                       end
                                                    end
                                                elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                    temp_E(1,m) = 1;
                                                    break
                                                end
                                            end
                                         end
                                    end
                                 end
                           end %end of build E and H
                           Z = [Z; temp_Z];
                           E = [E; temp_E];
                           H = [H; temp_H];
                           EQ = [EQ; zeros(1,h)];
                           F = [F; zeros(1,h)];
                           RHS = [RHS; temp_RHS];
                           if(k == 2 && plot_if_k_is_2)
                             colors = plot_RHS(temp_RHS,c,h,x,colors);
                           end
                           kappa = [kappa; l];
                           pivot_from_prev = [pivot_from_prev; 0,l,0];
                           previous_pivot = [previous_pivot; size(pivots_made,1)];
                           % build P and iota
                           if(basis_dim(row) == 1)
                                P = [P;row];
                                iota = [iota;l];
                           else
                                P = [P;P(row)];
                                iota = [iota;iota(row)];
                           end
                           % set up nlp_d
                           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                            x0 = feas_pt(1:k+1);
                           leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                           nonlcon = @(x) create_nonlcon(leq,[],x);
                           [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                           if(exitflag >= 1)
                             feas_pts = [feas_pts;sol];
                             basis_dim = [basis_dim;1];
                             bases = [bases;new_basis(1:h)];
                             regions = [regions;simplify(G_(:,end))];
                           elseif(exitflag == -2)
                             feas_pts = [feas_pts;x0];
                             basis_dim = [basis_dim;0];
                           else
                            feas_pts = [feas_pts;x0]; 
                            basis_dim = [basis_dim;0];
                            disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                            return  
                           end
                        end
                     else
                        for j = 1:h
                           if(j ~= l)
                               if(current_basis(j) == 0) % l and j are both w's
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [1,1];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j+h);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0];
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               else % l is a w, j is a z
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [1,0];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l+h),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l+h);
                                         G_ = pivot(G_,l,j);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+2),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0]; 
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end 
                               end
                           end
                        end
                     end
                  else % l is a z
                      if(G(l,l) ~= 0)
                        new_basis = current_basis;
                        new_basis(l) = 0;
                        if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0) 
                           new_basis(h+1) = size(Z,1)+1;
                           if(print_stuff)
                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                new_basis
                           end
                           basis_stack = [new_basis; basis_stack];
                           discovered_bases = [new_basis; discovered_bases];
                           G_ = pivot(G,l,l);
                           temp_Z = zeros(1,h);
                           temp_E = zeros(1,h);
                           temp_H = zeros(1,h^2);
                           temp_RHS = sym(zeros(1,h));
%                            new_pt = x0;
                           %build Z and RHS
                           for m = 1:h
                                [N,d] = numden(G_(m,end));
                                n = coeffs(N,'All');
                                if(isempty(n))
                                    temp_Z(1,m) = 1;
                                else
                                    val = subs(d,x,feas_pt);
                                    if(size(n,1) > 1 || abs(n(1)) > .00001)
                                        if(val > 0) 
                                            temp_RHS(1,m) = -N;
                                        else
                                            temp_RHS(1,m) = N;
                                        end
                                    end
                                end
                           end
                           %build E and H
                           for m = 1:h
                               one_was_feas = 0;
                                 if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                    for n = 1:h
                                        if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                            if(temp_H(1,(m-1)*h + n) ~= 1)
                                                fun = @(x) -x(k+1);
%                                                 x0 = feas_pt(1:k+1);
                                                temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                if(exitflag == 0)
                                                    disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                    return
                                                end
                                                if(exitflag >= 1 && -fval >= 0)
                                                    one_was_feas = 1;
                                                end
                                                if(exitflag >= 1 && abs(fval) < .00001)
                                                    temp_H(1,(m-1)*h + n) = 1;
                                                    for p = 1:h
                                                       if(temp_H(1,(n-1)*h + p) == 1) 
                                                           temp_H(1,(m-1)*h + p) = 1;
                                                       end
                                                    end
                                                elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                    temp_E(1,m) = 1;
                                                    break
                                                end
                                            end
                                         end
                                    end
                                 end
                           end %end of build E and H
                           Z = [Z; temp_Z];
                           E = [E; temp_E];
                           H = [H; temp_H];
                           EQ = [EQ; zeros(1,h)];
                           F = [F; zeros(1,h)];
                           RHS = [RHS; temp_RHS];
                           if(k == 2 && plot_if_k_is_2)
                             colors = plot_RHS(temp_RHS,c,h,x,colors);
                           end
                           kappa = [kappa; l];
                           pivot_from_prev = [pivot_from_prev; 0,l,0];
                           previous_pivot = [previous_pivot; size(pivots_made,1)];
                           % build P and iota
                           if(basis_dim(row) == 1)
                                P = [P;row];
                                iota = [iota;l];
                           else
                                P = [P;P(row)];
                                iota = [iota;iota(row)];
                           end
                           % set up nlp_d
                           fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                            x0 = feas_pt(1:k+1);
                           leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                           nonlcon = @(x) create_nonlcon(leq,[],x);
                           [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                           if(exitflag >= 1)
                             feas_pts = [feas_pts;sol];
                             basis_dim = [basis_dim;1];
                             bases = [bases;new_basis(1:h)];
                             regions = [regions;simplify(G_(:,end))];
                           elseif(exitflag == -2)
                             feas_pts = [feas_pts;x0];
                             basis_dim = [basis_dim;0];
                           else
                            feas_pts = [feas_pts;x0];   
                            basis_dim = [basis_dim;0];
                            disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                            return  
                           end
                        end
                     else
                        for j = 1:h
                           if(j ~= l)
                               if(current_basis(j) == 0) % l is a z and j is a w
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [0,1];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j+h),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l);
                                         G_ = pivot(G_,l,j+h);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(one_was_feas && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0]; 
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               else % l is a z, j is a z
                                  new_basis = current_basis;
                                  new_basis([l,j]) = [0,0];
                                  if( is_in(new_basis(1:h),discovered_bases(:,1:h)) == 0)
                                      if( subs(G(j,l),x,feas_pt) > 0 && subs(G(l,j),x,feas_pt) ~= 0)
                                         G_ = pivot(G,j,l);
                                         G_ = pivot(G_,l,j);
                                         G_([l,j],:) = G_([j,l],:);
                                         temp_Z = zeros(1,h);
                                         temp_E = zeros(1,h);
                                         temp_H = zeros(1,h^2);
                                         temp_RHS = sym(zeros(1,h));
%                                          new_pt = x0;
                                         %build Z and RHS
                                         for m = 1:h
                                            [N,d] = numden(G_(m,end));
                                            n = coeffs(N,'All');
                                            if(isempty(n))
                                                temp_Z(1,m) = 1;
                                            else
                                                val = subs(d,x,feas_pt);
                                                if(size(n,1) > 1 || abs(n(1)) > .00001)
                                                    if(val > 0) 
                                                        temp_RHS(1,m) = -N;
                                                    else
                                                        temp_RHS(1,m) = N;
                                                    end
                                                end
                                            end
                                         end
                                         %build E and H
                                         for m = 1:h
                                             one_was_feas = 0;
                                             if(temp_Z(1,m) ~= 1 && temp_E(1,m) ~= 1 && size(coeffs(temp_RHS(1,m),'All'),2) > 1)
                                                for n = 1:h
                                                    if(temp_Z(1,n) ~= 1 && temp_E(1,n) ~= 1 && n ~= m)
                                                        if(temp_H(1,(m-1)*h + n) ~= 1)
                                                            fun = @(x) -x(k+1);
%                                                             x0 = feas_pt(1:k+1);
                                                            temp = temp_RHS(1,setdiff(indices,[m,n]));
                                                            leq = matlabFunction([c,temp(temp~=0),temp_RHS(1,n)+x(k+1)],'Vars',x);
                                                            eq = matlabFunction(temp_RHS(1,m),'Vars',x);
                                                            nonlcon = @(x) create_nonlcon(leq,eq,x);
                                                            [~,fval,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,lb(1:k+1),ub(1:k+1),nonlcon,options,options2,options3);
                                                            if(exitflag == 0)
                                                                disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                                                return
                                                            end
                                                            if(exitflag >= 1 && -fval >= 0)
                                                                one_was_feas = 1;
                                                            end
                                                            if(exitflag >= 1 && abs(fval) < .00001)
                                                                temp_H(1,(m-1)*h + n) = 1;
                                                                for p = 1:h
                                                                   if(temp_H(1,(n-1)*h + p) == 1) 
                                                                       temp_H(1,(m-1)*h + p) = 1;
                                                                   end
                                                                end
                                                            elseif(one_was_feas == 0 && (exitflag == -2 || (exitflag >= 1 && -fval < 0)))
                                                                temp_E(1,m) = 1;
                                                                break
                                                            end
                                                        end
                                                     end
                                                end
                                             end
                                         end %end of build E and H
                                         %set up nlp_a
                                         discard1 = [];
                                         discard2 = [];
                                         for m = 1:h
                                            if(Z(row,m) == 1 || H(row,(l-1)*h + m) == 1 || m == l)
                                                discard1 = [discard1,m];
                                            end
                                            if(temp_Z(1,m) == 1 || temp_H(1,(j-1)*h + m) == 1 || m == j)
                                                discard2 = [discard2,m];
                                            end
                                         end
                                         fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                          x0 = feas_pt(1:k+1);
                                         temp = [RHS(row,setdiff(indices,discard1)), temp_RHS(1,setdiff(indices,discard2))]+x(k+1);
                                         leq = matlabFunction([c+x(k+1),temp(temp~=0)],'Vars',x);
                                         eq = matlabFunction(RHS(row,l),'Vars',x);
                                         nonlcon = @(x) create_nonlcon(leq,eq,x);
                                         [~,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                         if(exitflag == 0)
                                            disp(sprintf('Fmincon failed. Line %d',MFileLineNr()));
                                            return
                                         elseif(exitflag >= 1)
                                             new_basis(h+1) = size(Z,1)+1;
                                             if(print_stuff)
                                                disp(sprintf('Adding new basis. Line %d',MFileLineNr()));
                                                new_basis
                                             end
                                             basis_stack = [new_basis; basis_stack];
                                             discovered_bases = [new_basis; discovered_bases];
                                             Z = [Z; temp_Z];
                                             E = [E; temp_E];
                                             H = [H; temp_H];
                                             EQ = [EQ; zeros(1,h)];
                                             F = [F; zeros(1,h)];
                                             RHS = [RHS; temp_RHS];
                                             if(k == 2 && plot_if_k_is_2)
                                                colors = plot_RHS(temp_RHS,c,h,x,colors);
                                             end
                                             kappa = [kappa; j];
                                             pivot_from_prev = [pivot_from_prev; 1,l,j];
                                             previous_pivot = [previous_pivot; size(pivots_made,1)];
                                             % build P and iota
                                             if(basis_dim(row) == 1)
                                                P = [P;row];
                                                iota = [iota;l];
                                             else
                                                P = [P;P(row)];
                                                iota = [iota;iota(row)];
                                             end
                                             % set up nlp_d
                                             fun = @(x) 0; if(use_optimal_instead_of_feasible)     fun = @(x) -x(k+1); end
%                                              x0 = feas_pt(1:k+1);
                                             leq = matlabFunction([c+x(k+1),temp_RHS(temp_RHS~=0) + x(k+1)],'Vars',x);
                                             nonlcon = @(x) create_nonlcon(leq,[],x);
                                             [sol,~,exitflag] = use_fmincon(fun,x0,A,b,Aeq,beq,[lb(1:k),.0001],ub(1:k+1),nonlcon,options,options2,options3);
                                             if(exitflag >= 1)
                                                 feas_pts = [feas_pts;sol];
                                                 basis_dim = [basis_dim;1];
                                                 bases = [bases;new_basis(1:h)];
                                                 regions = [regions;simplify(G_(:,end))];
                                             elseif(exitflag == -2)
                                                 feas_pts = [feas_pts;x0];
                                                 basis_dim = [basis_dim;0];
                                             else
                                                feas_pts = [feas_pts;x0]; 
                                                basis_dim = [basis_dim;0];
                                                disp(sprintf('Fmincon failed. Exitflag: %d, Line %d',exitflag,MFileLineNr()));
                                                return  
                                             end
                                         end                                     
                                      end
                                  end
                               end
                           end
                        end
                     end
                  end
              end
           end %end find adjacent
           %eliminate overlapping adjacent
%            for i = 1:h
%               if(Z(row,i) == 1)
%                  if(current_basis(i) == 0)
%                      if(G(i,i+h) ~= 0)
%                         new_basis = current_basis;
%                         new_basis(i) = 1;
%                         discovered_bases = [new_basis; discovered_bases];
%                      end
%                  else
%                      if(G(i,i) ~= 0)
%                         new_basis = current_basis;
%                         new_basis(i) = 0;
%                         discovered_bases = [new_basis; discovered_bases];
%                      end
%                  end
%               end
%            end
       end
   end
   previous_basis = current_basis;
end

% disp(sprintf('Phase 2 iteration: %d',phase2_iter));

fileID = fopen('solution.txt','w');
for i = 1:size(bases,1)
    fprintf(fileID,'Region %d:\n\n',i);
    for j = 1:h
        if(bases(i,j) == 0)
            fprintf(fileID,'w%d:  ',j);
        else
            fprintf(fileID,'z%d:  ',j);
        end
        fprintf(fileID,'%s\n',regions(h*(i-1)+j));
    end
    fprintf(fileID,'\n\n\n');
end
fclose(fileID);
warning('on','MATLAB:rankDeficientMatrix');