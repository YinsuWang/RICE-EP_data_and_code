function w = welfareMugrowth(argument,s,Para,Var,Initial)
    mu = growth2013(argument(1,:),argument(2,:),argument(3,:),Para.Tmax);
    mu(mu>1) = 1;
    mu(mu<0) = 0;
    result = OBJMuNWU(mu,s,Para,Var,Initial);
    w = sum(result.wel,"all"); 
end