function [c,ceq] = tempconMuNWU(argument,s,Para,Var,Initial)
    [~,result] = welfare_Mu_NWU(argument,s,Para,Var,Initial);

 c = result.T(size(result.T),1)-2;
 %c = max((temp.T),[],"all")-2;
 ceq =[];
end