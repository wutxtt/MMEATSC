function [Population,Species,center] = Species_division(Population,Nns,mod, Global)
% The environmental selection of NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% normalization in the decision space 对决策空间的标准化
    z=min(Population.decs,[],1);
    Z=max(Population.decs,[],1);
    
    decc=(Population.decs-z)./repmat(Z-z,length(Population),1);
    %% 决策空间自适应半径
    NP = length(Population); %种群规模
    NB = min(Nns,max(NP-1,0)); %取有最小距离个体的个数
%     decc = Population.decs;
    d_dec = pdist2(decc,decc,'euclidean');  %计算决策空间每个点之间的距离
    temp = sort(d_dec); %对于Nns个个体之间的距离进行按列升序排序
    r_dec = sum(sum(temp(1:NB+1,:)))./(NB.*NP); %计算NB*NP个距离的平均距离
    [n,dim]=size(Population.decs);
    S = (1:n);  % S为1*n的矩阵
    k = 1;
    center = [];
    while size(S) > 0
        s = randperm(length(S));
        x = S(s(1));  % 随机选择一个未被划分物种的个体作为物种的中心个体 
        % C为包含全部物种的集合
        C(k).p = [];  % 初始化一个物种中个体的list
        Distance = d_dec(x,:);  % 取这个中心个体与全部其他个体在决策空间的距离
        Nbs = find(Distance < r_dec);  % 得到所有在决策空间距离小于r_dec的个体
        C(k).p = [C(k).p Nbs];  % 将所得到的个体全部放入初始化的list中
        % 遍历C中之前所有的物种
        i = k - 1; 
        while i > 0
            C(k).p = setdiff(C(k).p,C(i).p);  % 对当前物种中的个体对所有之前的物种进行去重，返回在k中存在但不在i中存在的元素
            i = i - 1;
        end
        S = setdiff(S,C(k).p);  % 将已经被划分进入物种中的个体从原始种群中分离，返回S中存在但不在k中存在的元素
        center = [center x];  % 将本次循环的物种中心个体放入center中
        k = k + 1; % 物种个数+1
    end
    center = Population(center);
    %% normalization in the decision space 对决策空间的标准化
    z=min(Population.objs,[],1);
    Z=max(Population.objs,[],1);
    
    objj=(Population.objs-z)./repmat(Z-z,length(Population),1);
    
    
    AllSpecies = [];
%     Seed = Population(Seed);
    if mod == 1
         for i = 1:length(C)
%             beta = 5;
%                 if length(Species(i).p) >= beta
                    Species(i).p = Population(C(i).p);
%                 end
        end
        Population = Population;
    end
    if mod == 2
        size_population = 0;
        for i = 1:length(C)
            Species(i).p = Population(C(i).p);
            size_population = size_population + length(Species(i).p);
        end
        
        SeedSet = [];
        for i = 1:length(Species)
%             for j = 1:length(Species(i).p)
            nor_obj = (Species(i).p.objs-min(Species(i).p.objs))./(max(Species(i).p.objs)-min(Species(i).p.objs));
            sum_obj = sum(nor_obj,2);
            [~,R] = sort(sum_obj);
            SeedSet = [SeedSet Species(i).p(R(1))];
        end
        [FrontNo,MaxFNo] = NDSort(SeedSet.objs,inf);
        tabulate(FrontNo(:))
        
        
        while size_population > Global.N
            num = [];
            for i = 1:length(Species)
                num = [num length(Species(i).p)];
            end
             m = find(num == max(num));
             s = randperm(length(m));
             x = m(s(1));
             
             [FrontNo,MaxFNo] = NDSort(Species(x).p.objs,inf);
             Last     = find(FrontNo==MaxFNo); %最后非支配集的位置索引
             CrowdDis = Crowding(Species(x).p(Last).objs);
%              CrowdDis = Crowding(Species(x).p.objs);
             [~,R] = sort(CrowdDis);
             Species(x).p(Last(R(1))) = []; %最后非支配层剔除位置索引    
%              Species(x).p(R(1)) = []; %最后非支配层剔除位置索引    

%              CrowdDis = Crowding(Species(x).p.objs);
% %              CrowdDis = Crowding(Species(x).p.objs);
%              [~,R] = sort(CrowdDis);
%              Species(x).p(R(1)) = []; %最后非支配层剔除位置索引    
             size_population = size_population - 1;
        end
        Population = [];
        for i = 1:length(Species)
            Population = [Population Species(i).p];
        end
    end
end