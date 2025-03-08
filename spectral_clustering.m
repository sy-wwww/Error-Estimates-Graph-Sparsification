function er1B =spectral_clustering(n, cedge, B, r)
% Algorithm 2 for spectral clustering from the paper "Empirical Error Estimates for Graph Sparsification"
% Inputs:
%   - n: Number of nodes in the graph
%   - B: Bootstrapping sample size for the first loop
%   - cedge: An N×4 matrix where:
%       * The first two columns represent the nodes of the sampled edges.
%       * The third column contains the weights of the edges in the sparsified graph.
%       * The fourth column indicates the number of times each edge is sampled.
%   - r: A number that the user believes is safely above the correct number of clusters

    S0 = sparse(cedge(:,1),cedge(:,2),cedge(:,3),n,n);
    S = S0 + transpose(S0);
    SS = sum(S);
    S = diag(-SS) + S;
    k = sum(cedge(:,4));
    k1 = length(cedge(:,1));


    [uSM,d,v] = eigs(full(S),r,"smallestreal");


    eigv = flip(diag(d));
    eigv = eigv(1:(r-1));
    

    er1B = zeros(B,1);

    for b = 1:B
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


        [uS1M,d,v] = eigs(full(S1),r,"smallestreal");

        eigvb = flip(diag(d));
        eigvb = eigvb(1:(r-1));


        er1B(b) = max(abs(eigvb./eigv - 1));
    end

end
    
   