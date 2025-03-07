function f = func(L,S,lambda,y)
    f1 = norm(S-L,"fro");
    f2 = svds(S-L,1);
    n = length(y);
    betaL = (diag(ones(n,1)) + lambda*L)\y;
    betaS = (diag(ones(n,1)) + lambda*L)\y;
    f3 = norm(betaS-betaL);
    f = [f1,f2,f3];
end
