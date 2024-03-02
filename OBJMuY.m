function result = OBJMuY(mu,s,Para,Var,Initial)

Tmax = Para.Tmax;
I = Para.I;
N = Para.N;     %前面三个在main中确定
tau = Para.tau/100*N; %碳税收系数
d = Para.d/100*N;     %气候损害分布
q = Para.q/100*N;     %消费分布，每个收入组的消费占比
ryc = Para.ryc;
HT = Para.HT;
CT = Para.CT;
%Tem = Initial.Tem;
%Regions = Para.Regions;
%Var.Tsd = sqrt(365.25)*Var.Tsd;
%Transfer exogenous variables
L = Var.L;              %时间*收入组*区域
Ltol = permute(sum(L,2),[1,3,2]);        %时间*区域 2是维度
cbtol = permute(sum(Var.cb,2),[1,3,2]);
A = Var.A;              %技术水平
sigma = Var.sigma;      %碳强度
Eland = Var.EL0;      %非工业碳排放

%只有消费Consumption、福利Walfare是三维，资本Capital等还是二维
    %Consumption
    c = zeros(Tmax,N,I);
    %control variable
    %mu = X(:,I+1:2*I) ; %将mu同样作为控制变量 emissions reduction rate
    %mu = mu.*ones(Tmax,I);
    %s = X(:,1:I);       %saving rate
    s = s.*ones(Tmax,I);
    S = zeros(Tmax,I);     %Total savings
    %考虑下capital的问题，是否是控制变量，在RCK model中social planner's problem 一般都是 allocate k c，
    % 但是DICE/RICE中不一定，许多模型（尤其是考虑碳税）都将储蓄率设为常数
    wel = zeros(Tmax,N,I);
    %第一期的Q可以考虑直接用GDP，k(1,:)为上一年投资存量
    K(1,:) = Initial.K;
    Tax = mu.^(Para.theta(2,1)-1).*Var.pb/1000;% trillion US$ per GtCO2
    
    %预设部分空数组
    E = zeros(Tmax,1);
    F = E;
    T = E;
    Tlo = E;
    %Var.O = Var.O*0;
    %Eland = Eland*0;
    
    AT = 1;% atmosphere,
    UP = 2;% upper oceans
    LO = 3;% biosphere
    CO2 = zeros(Tmax,3);
    
    Q(1,:) = A(1,:).*K(1,:).^(Para.alpha).*Ltol(1,:).^(1-Para.alpha);%注意period，L(1)、A(1)对应的是第0期,已经在ExogenousVariables中处理
    S(1,:) = Q(1,:).*s(1,:);
    K(2,:) = (1-Para.delta_K)*K(1,:) + S(1,:);%.*Ltol(1,:);%这里s可以作为控制变量(这里的S已经是总储蓄了，因为前面的Q已经是总GDP，所以不该×L）
    Eind(1,:) = sigma(1,:).*(1-mu(1,:)).*Q(1,:);
    Etol(1,:) = Eind(1,:) ;
    E(1) = sum(Etol(1,:),2)+ Eland(1);%CO2 emissions from land use, land use change and forestry.tCO2
    %E(1)=0;
    CO2(1,AT) = E(1)*(1/3.666)+(Para.phi(1,1)).*Initial.CO2(1,AT)+(Para.phi(2,1)).*Initial.CO2(1,UP);%DICE 中为5年，所以是5/3.66
    CO2(1,UP) = (Para.phi(1,2)).*Initial.CO2(1,AT)+(Para.phi(2,2)).*Initial.CO2(1,UP)+(Para.phi(3,2)).*Initial.CO2(1,LO);
    CO2(1,LO) = (Para.phi(2,3)).*Initial.CO2(1,UP)+(Para.phi(3,3)).*Initial.CO2(1,LO);
    F(1) = Para.eta*(log(CO2(1,AT)/Initial.CO2pre(1,1))/log(2))+Var.O(1);
    %F(1) = 0;
    T(1) = Initial.T+Para.xi(1)*(F(1)-Para.xi(2)*Initial.T-Para.xi(3)*(Initial.T-Initial.Tlo));
    Tlo(1) = Initial.Tlo + Para.xi(4)*(Initial.T-Initial.Tlo);

    D(1,:) = (T(1)+Initial.T)*Para.theta(3,:)+(T(1)+Initial.T)^2*Para.theta(4,:);
    Lambda(1,:) = Var.theta1(1,:).*mu(1,:).^Para.theta(2,1);%这里可能有矩阵大小问题
    
    Y(1,:) = (1-Lambda(1,:))./(1.+D(1,:)).*Q(1,:);
    

    for t = 2:Tmax
        
        Q(t,:) = A(t,:).*K(t,:).^(Para.alpha).*Ltol(t,:).^(1-Para.alpha);%注意period，L(1)、A(1)对应的是第0期,已经在ExogenousVariables中处理
        S(t,:) = Q(t,:).*s(t,:);
        K(t+1,:) = (1-Para.delta_K)*K(t,:) + S(t,:);%.*Ltol(t,:);
        %delmu =  mu(t,:)-mu(t-1,:);
        %mu(t,(delmu<0)) = mu(t-1,(delmu<0));
        Eind(t,:) = sigma(1,:).*(1-mu(t,:)).*Q(t,:);
        Etol(t,:) = Eind(t,:)  ;%+ Eres(t,:);
        E(t) = sum(Etol(t,:),2)+ Eland(t);
        %E(t) = 0;
        CO2(t,AT) = E(t)*(1/3.666)+(Para.phi(1,1)).*CO2(t-1,AT)+(Para.phi(2,1)).*CO2(t-1,UP);%DICE 中为5年，所以是5/3.66
        CO2(t,UP) = (Para.phi(1,2)).*CO2(t-1,AT)+(Para.phi(2,2)).*CO2(t-1,UP)+(Para.phi(3,2)).*CO2(t-1,LO);
        CO2(t,LO) = (Para.phi(2,3)).*CO2(t-1,UP)+(Para.phi(3,3)).*CO2(t-1,LO);
        F(t) = Para.eta*(log(CO2(t,AT)/Initial.CO2pre(1,1))/log(2))+Var.O(t);
        %F(t) = 0;
        T(t) = T(t-1)+Para.xi(1)*(F(t)-Para.xi(2)*T(t-1)-Para.xi(3)*(T(t-1)-Tlo(t-1)));
        Tlo(t) = Tlo(t-1) + Para.xi(4)*(T(t-1)-Tlo(t-1));

        D(t,:) = (T(t)+Initial.T)*Para.theta(3,:)+(T(t)+Initial.T)^2*Para.theta(4,:);
        Lambda(t,:) = Var.theta1(t,:).*mu(t,:).^Para.theta(2,1);
        
        Y(t,:) = (1-Lambda(t,:))./(1.+D(t,:)).*Q(t,:);
    end
    K(size(K,1),:) = [];
    C = (Y-S)./Ltol;
    %[HDD,CDD] = DegreeDaysThoms(HT,CT,T +Initial.Tpre + Var.Tdel,Var.Tsd,365.25);
    [HDD,CDD] = DegreeDaysGelegenis(HT,CT,T +Initial.Tpre+ Var.Tdel,Var.Temmax,Var.Temmin,365.25);
    HDD(HDD<0) = 0;
    CDD(CDD<0) = 0;
    %Tem(:,3) = array2table(table2array(Tem(:,3))+(T(1)-Initial.T));
    %[HDD(1,:),CDD(1,:)] = DegreeDays(Tem,Regions,'Temperature','Region','date',HT,CT) ;
    %Tsd:各地区温度的标准差，暂时设定不随时间变化
        %Tdel 也需要通过Region.m来处理
    ECH = Para.gamma(1,1).*HDD + Para.gamma(1,2);%Energy consumption percapital of heating
    ECC = Para.gamma(2,1).*CDD + Para.gamma(2,2);%Energy consumption percapital of cooling
    %ECH = exp(Para.gamma(1,1).*log(HDD) + Para.gamma(1,2));
    %ECC = exp(Para.gamma(1,1).*log(CDD) + Para.gamma(1,2));
    Eres = ECH + ECC;%人均residential 能源消费

    EXPH = Para.price(1,:).*ECH; %取暖能源可以包括煤、天然气、电力，这里采用居民全部能源能源消费的单位价格
    EXPC = Para.price(2,:).*ECC; %可以认为制冷能源仅有电力，这里直接使用电力价格，
    EXPH(EXPH<0) = 0;
    EXPC(EXPC<0) = 0;
    ce = EXPH+EXPC;
    cpre = (Q - S)./Ltol - ce ; %可以假设cb（基本消费）为0,cbar 应为人均消费
    %cbar = (Y - S)./Ltol - ce ;
    cbar = (Q - S)./Ltol;
    % = bsxfun(@times,cpre,q)-bsxfun(@times,cbar.*D,d)-bsxfun(@times,(cpre.*Lambda+Tax.*Etol./Ltol),tau)+bsxfun(@times,Tax.*Etol./Ltol,ryc);%考虑碳税后的消费
    for t = 1:Tmax
        GrossConsumption = bsxfun(@times,cpre(t,:),q);
        DamageCost = bsxfun(@times,cbar(t,:).*D(t,:),d);
        MitigationCost = bsxfun(@times,(cpre(t,:).*Lambda(t,:)),tau);
        TaxPayment = bsxfun(@times,(Tax(t,:).*Etol(t,:)./Ltol(t,:)),tau);
        Refund = bsxfun(@times,Tax(t,:).*Etol(t,:)./Ltol(t,:),ryc);
        ctol(t,:,:) = GrossConsumption - DamageCost - MitigationCost - TaxPayment + Refund;
        %ctol(t,:,:) = bsxfun(@times,cpre(t,:),q)-bsxfun(@times,cbar(t,:).*D(t,:),d)-bsxfun(@times,(cpre(t,:).*Lambda(t,:)),tau)-bsxfun(@times,(Tax(t,:).*Etol(t,:)./Ltol(t,:)),tau)+bsxfun(@times,Tax(t,:).*Etol(t,:)./Ltol(t,:),ryc);%考虑碳税后的消费
        Ep(t,:,:) = ce(t,:)./bsxfun(@times,Y(t,:),q);
        Ep(t,:,:) = Ep(t,:,:).*L(t,:,:);
    end
    CE = zeros(Tmax,N,I);%生成一个三维的数组用于后面的计算
    for i = 1:N
        CE(:,i,:) = ce;
    end
    c = ctol - CE;%将能源消费视为基本支出的超额消费
    %Ep(1,:,:) = (1-c(1,:,:)-Var.cb(1,:,:))./(ce+Var.cb(1,:,:)+cpre);
    Eptol = ce./(Y./Ltol);
    x=-linspace(1,Para.Tmax,Para.Tmax);
    R = ((1+Para.beta).^(x))';
    if Para.epsilon == 1
       U = log(c).*L;
    else
       U = c.^(1-Para.epsilon)./(1-Para.epsilon);
    end
    wel = times(R,ones(1,5,12)).*U;

    Ctrue = R*ones(1,12).*(C - ce)./(Y./Ltol);
    %Saving results
    %可能需要调整,在求最优解的时候，function只能由一个输出值，目标函数为c/W，因此可能无法通过本程序输出结果
    varlist = who;
    for i = 1:size(varlist,1)
        result.(varlist{i}) = eval([(varlist{i})]);
    end

    
