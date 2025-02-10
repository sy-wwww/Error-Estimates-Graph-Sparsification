function f = func(L,S)
    f1 = norm(S-L,"fro");
    f2 = svds(S-L,1);
    f = [f1,f2];
end