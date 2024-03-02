function [w,result] = welfare_Mu(argument,s,Para,Var,Initial)

if size(argument,1) ==1
    if size(argument,2) ==1
        mu = argument*ones(Para.Tmax,Para.I);
    else
        mu = ones(Para.Tmax,1)*argument;
    end
elseif size(argument,1)<Para.Tmax
    if size(argument,2) ==1
    mu = [argument;ones(Para.Tmax-size(argument,1),1)]*ones(1,Para.I);
    elseif size(argument,2) ==12
    %mu = zeros(Para.Tmax,Para.I);
    %Tax0 = argument;
    [g,delta ]=decline2013(argument,Para.gAmin,Para.deltaAmin);
    mu0= growth2013(argument(size(argument,1),:),g,delta,Para.Tmax-size(argument,1));
    mu0(1,:) =[];
    %mu0 = ones(Para.Tmax-size(argument,1),Para.I);
    mu = [argument;mu0];
    %mu(1:size(Tax0,1),1:size(Tax0,2)) = Tax0;
    end
else
    mu = argument;
end

    result = OBJMuY(mu,s,Para,Var,Initial);
    w = sum(result.c,"all"); 
end