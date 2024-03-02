function [c,ceq] = tempconMU(argument)
    global gs gPara gVar gInitial  ;
    [~,result] = welfare_Mu(argument,gs,gPara,gVar,gInitial);

 c = result.T(size(result.T),1)-2;
 %c = max((temp.T),[],"all")-2;
 ceq =[];
end