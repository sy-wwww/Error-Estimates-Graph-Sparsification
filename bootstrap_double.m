function [er1B, er1B_sd] = bootstrap_double(n, cedge, B, C, num_func)
% Algorithm 1 from the paper "Empirical Error Estimates for Graph Sparsification"
% Inputs:
%   - n: Number of nodes in the graph
%   - B: Bootstrapping sample size for the first loop
%   - C: Bootstrapping sample size for the second loop
%   - num_func: Number of functions to be estimated
%   - cedge: An N×4 matrix where:
%       * The first two columns represent the nodes of the sampled edges.
%       * The third column contains the weights of the edges in the sparsified graph.
%       * The fourth column indicates the number of times each edge is sampled.

%     cedge = zeros(k1,3);
%     cedge(:,1:2) = nonzero1(a(:,1),:);
%     cedge(:,3) = -L0(nonzero(a(:,1)))./p1(a(:,1))/k.*a(:,2);
    S0 = sparse(cedge(:,1),cedge(:,2),cedge(:,3),n,n);
    S = S0 + transpose(S0);
    SS = sum(S);
    S = diag(-SS) + S;
    k = sum(cedge(:,4));
    k1 = length(cedge(:,1));
    
    er1B = zeros(B,num_func);
    er1B_sd = zeros(B,num_func);
    parfor b = 1:B
        s1 = randsample(1:k1,k,true,cedge(:,4)/k);
        a1 = tabulate(s1);
        a1 = transpose(a1(a1(:,2)>0,1:2));
        k2 = length(a1(1,:));

        cedge_b = zeros(k2,3);

        cedge_b(:,1:2) = cedge(a1(1,:),1:2);
        cedge_b(:,3) = cedge(a1(1,:),3)./cedge(a1(1,:),4) .* transpose(a1(2,:));

        S0 = sparse(cedge_b(:,1),cedge_b(:,2),cedge_b(:,3),n,n);
        S1 = S0 + transpose(S0);
        SS = sum(S1);
        S1 = diag(-SS) + S1;
        %%%%%%%%%%%%%%%%%%%%%

        er1B1 = zeros(C,num_func);
        for t = 1:C
            s2 = randsample(1:k2,k,true,a1(2,:)/k);
            a2 = tabulate(s2);
            a2 = transpose(a2(a2(:,2)>0,1:2));
            k3 = length(a2(1,:));

            cedge_b1 = zeros(k3,3);
            cedge_b1(:,1:2) = cedge_b(a2(1,:),1:2);
            cedge_b1(:,3) = cedge_b(a2(1,:),3)./transpose(a1(2,a2(1,:))) .* transpose(a2(2,:));
            S01 = sparse(cedge_b1(:,1),cedge_b1(:,2),cedge_b1(:,3),n,n);
            S11 = S01 + transpose(S01);
            SS1 = sum(S11);
            S11 = diag(-SS1) + S11;

            er1B1(t,:) = func(S1,S11);

        end

        val = func(S,S1);
        
        er1B(b,:) = val;
        er1B_sd(b,:) = (val - mean(er1B1))./std(er1B1).*std(er1B)+mean(er1B);
    end
end

