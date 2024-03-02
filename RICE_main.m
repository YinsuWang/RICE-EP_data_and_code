%RICE_EPmain
tic
clear
warning("off");
for Negishi=0:1
    %初始参数设置
    %SetP.Step=5;%设定每几期做一次决定
    SetP.Negishi = Negishi; %i=0 No Negishi welfare weight;
    SetP.Tend = 2100;
    SetP.Tstart = 2015;%考虑设置为2015
    SetP.HT = 18.3;  %设定平衡点温度用于计算HDD/CDD
    SetP.CT = 23.9;
    SetP.Tmax = SetP.Tend-SetP.Tstart;%Tmax需要根据后续数据重新确认
    %SetP.I = size(unique(CR(:,2)),1);%根据数据确定数组大小 
    SetP.N = 5;%设置收入组数
    %设置A L 的校准参数
    SetP.Tpop = 0;  %采用近Tpop年的数据进行校准，Tpop=0则使用外部数据（WPP）
    SetP.TA = 20;    %采用近TA年的数据进行校准
    
    try
        load("InputData.mat");
    catch ERROR
    end
    
    %设置碳税返还场景
    scenarios = ["NoTax","NoRecycling","RecyclingTax","RecyclingCost"];
    goals = ["","2C"];
    nonlcon = [];
    Algorithm = 'interior-point';
    init = 0.2*ones(ceil(Para.Tmax) ,Para.I);
    for k = 1:2
    for n = 1:size(scenarios,2)
    scenario = scenarios(1,n);
    goal = goals(k);
    if scenario == "NoTax"
        Para.tau = ones(Para.N,Para.I)*20;
        Para.ryc = Para.tau;%此处tau为0.2，对应的是总消费，若在人均消费/损害等使用，应为1
    elseif scenario == "NoRecycling"
        Para.ryc = zeros(Para.N,Para.I);
        Para.tau = Para.q;
    elseif scenario == "RecyclingTax"
        Para.ryc = ones(Para.N,Para.I)*20;
        Para.tau = Para.q;
    elseif scenario == "RecyclingCost"
        Para.ryc = ones(Para.N,Para.I)*20;
        Para.tau = zeros(Para.N,Para.I);
    elseif scenario == "DC2"
        nonlcon = @(argument)tempconMuNWU(argument,s,Para,Var,Initial);
        Para.tau = Para.q;
        init = optiargument;
    end
        if k ==1
            nonlcon = [];   
            if n>1
                init =  optiargument;
            end
        elseif k==2
            nonlcon = @(argument)tempconMuNWU(argument,s,Para,Var,Initial);
            init = Aresults.(scenarios(n)).mu(Para.period-Para.Tstart,:);
        end
    scenarioname = strcat(scenario,goal);
    
        Para.beta  = 0.015;
        Para.epsilon = 1;
          
        InitialValue.s = 0.258;
        options = optimset('Display','iter-detailed','MaxFunEvals',1000000,'Algorithm',Algorithm,'TolX',1e-8,'TolFun',1e-12,'MaxIter',10000,'UseParallel',true); %various other options
    
    s = InitialValue.s*ones(Para.Tmax,Para.I);
    
    %设定控制变量mu的形式以及初始值
    argument = init;
    argumentMU = zeros(3,Para.I);
      obj = @(argument)(-welfare(argument,s,Tax,Para,Var,Initial));
      nonlconMU = @tempconMU;
      objMU = @(argument)(-welfare_Mu(argument,s,Para,Var,Initial));
      ObjMuNWU = @(argument)(-welfare_Mu_NWU(argument,s,Para,Var,Initial));
      nonlconMuNWU = @tempconMuNWU;

    [optiargument,W,exitflag] = fmincon(ObjMuNWU,init,[],[],[],[],zeros(size(init)),ones(size(init)),nonlcon,options);

    optimu=optiargument;
    if k==2
    inits.(scenario) = optiargument;
    end

    [w,Aresults.(scenarioname)] = welfare_Mu_NWU(optimu,s,Para,Var,Initial);
    Aresults.(scenarioname).exitflag = exitflag;
    end
    end
    [w,Aresults.NoMitigation] =welfare_Mu_NWU(0,s,Para,Var,Initial);
    filename=['result',datestr(now,30),'.mat'];
    elapsedTime = toc;
    save(filename,'Aresults','elapsedTime');

    toc
end


