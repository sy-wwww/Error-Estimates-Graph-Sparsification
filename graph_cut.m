function er1B = graph_cut(n, cedge, B, CCS)
% Algorithm 2 from the paper "Empirical Error Estimates for Graph Sparsification"
% Inputs:
%   - n: Number of nodes in the graph
%   - B: Bootstrapping sample size
%   - cedge: An N×4 matrix where:
%       * The first two columns represent the nodes of the sampled edges.
%       * The third column contains the weights of the edges in the sparsified graph.
%       * The fourth column indicates the number of times each edge is sampled.
%   - CCS: A sparse matrix where each row is an n-dimensional binary vector representing a cut in the graph.

    S0 = sparse(cedge(:,1),cedge(:,2),cedge(:,3),n,n);
    S = S0 + transpose(S0);
    SS = sum(S);
    S = diag(-SS) + S;
    k = sum(cedge(:,4));
    k1 = length(cedge(:,1));

    CSD = diag(CCS * S *transpose(CCS));
    
    SD2 = CCS(:,nonzero1(s,1))+ CCS(:,nonzero1(s,2));
    SD2(find(SD2==2))=0;
    R = L0(nonzero(s))./p1(s);
    estd = std(transpose(SD2 .* transpose(R)))/sqrt(k);


    er1B = zeros(B,10);
    parfor b = 1:B
        s1 = randsample(1:k1,k,true,cedge(:,4)/k);
        a1 = tabulate(s1);
        a1 = transpose(a1(a1(:,2)>0,1:2));
        k2 = length(a1(1,:));

        cedge_b = zeros(k2,3);

        cedge_b(:,1:2) = cedge(a1(1,:),1:2);
        cedge_b(:,3) = cedge(a1(1,:),3)./a(a1(1,:),2) .* transpose(a1(2,:));

        S0 = sparse(cedge_b(:,1),cedge_b(:,2),cedge_b(:,3),n,n);
        S1 = S0 + transpose(S0);
        SS = sum(S1);
        S1 = diag(-SS) + S1;

        CS1D = diag(CCS * S1 *transpose(CCS));

        er1B(b) = max(abs(CSD-CS1D)./transpose(estd));
    end

end
