function graph_cut(N, cedge, B, C, CCS)
    S0 = sparse(cedge(:,1),cedge(:,2),cedge(:,3),N,N);
    S = S0 + transpose(S0);
    SS = sum(S);
    S = diag(-SS) + S;
    k = sum(cedge(:,4));
    k1 = length(cedge(:,1));
    
    CS = CCS * S *transpose(CCS);
    CSD = diag(CS);
    
    SD2 = CCS(:,nonzero1(s,1))+ CCS(:,nonzero1(s,2));
    SD2(find(SD2==2))=0;
    R = L0(nonzero(s))./p1(s);
    SD3 = SD2 .* transpose(R); 
    estd = std(transpose(SD3))/sqrt(k);


    er1B = zeros(B,10);
    parfor b = 1:B
        % [nonzeroS,nonzero1S, p1S] = sampling_probability_unweighted(S, -S0 - transpose(S0), -S0, method1);
        s1 = randsample(1:k1,k,true,cedge(:,4)/k);
        a1 = tabulate(s1);
        a1 = transpose(a1(a1(:,2)>0,1:2));
        k2 = length(a1(1,:));

        cedge_b = zeros(k2,3);

        cedge_b(:,1:2) = cedge(a1(1,:),1:2);
        cedge_b(:,3) = cedge(a1(1,:),3)./a(a1(1,:),2) .* transpose(a1(2,:));

        S0 = sparse(cedge_b(:,1),cedge_b(:,2),cedge_b(:,3),N,N);
        S1 = S0 + transpose(S0);
        SS = sum(S1);
        S1 = diag(-SS) + S1;
        %%%%%%%%%%%%%%%%%%%%%%%%

        % SD1 = CCS(:,nonzero1(s1,1))+ CCS(:,nonzero1(s1,2));
        % SD2 = SD1;
        % SD2(find(SD2==2))=0;
        % % R = L0(nonzero(s))./p1(s);
        % SD3 = SD2 .* transpose(R); 
        % estd = std(transpose(SD3))/sqrt(k);

        CS1 = CCS * S1 *transpose(CCS);
        CS1D = diag(CS1);
        CS1DM = min(CS1D);

        es1 = CCS(1,:) * S1 *transpose(CCS(1,:));


        er1B(b) = max(abs(CSD-CS1D)./transpose(estd));

    end
    boot_calc_stat(:,:,i) = get_stat(er1B);
    % er2_sd(i,:) = (er2(i,:)-mean(er1B(:,4:10)))./std(er1B(:,4:10));
end

    
    
    boot_calc_avg = mean(boot_calc_stat,3);
    boot_calc_sd = std(boot_calc_stat,0,3);
    
    
    
    % pop_calc1 = get_stat(er2);
    % pop_calc2 = get_stat(er2_sd);
    pop_calc2 = get_stat(er2);
    
    

end