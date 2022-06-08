function [Population, Seed, Species, Front] = Seed_detection(Population,Species,N,mod, Global,center)
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
    z=min(Population.objs,[],1);
    Z=max(Population.objs,[],1);
    objj=(Population.objs-z)./repmat(Z-z,length(Population),1);
    %% normalization in the decision space 对决策空间的标准化
    z=min(Population.decs,[],1);
    Z=max(Population.decs,[],1);
    decc=(Population.decs-z)./repmat(Z-z,length(Population),1);
    
    
    Seed = [];
    Sp_size = length(Population);
    min_Sp = 1;
    min_Sp_size = N;
    max_Sp = 1;
    max_Sp_size = 0;
    beta = 5;
    for i = 1:length(Species)
        if length(Species(i).p) < min_Sp_size
            min_Sp = i;
            min_Sp_size = length(Species(i).p);
        end
        if length(Species(i).p) > max_Sp_size
            max_Sp = i;
            max_Sp_size = length(Species(i).p);
        end
            [FrontNo,MaxFNo] = NDSort(Species(i).p,inf);
            First = find(FrontNo==1);
            Seed(i).p = [];
            Seed(i).p = Species(i).p(First);
    end
    
    if mod == 1
        Seed_temp = [];
        Front = [];
        for i = 1:length(Seed)
            NewSeedSet(i).p = Seed(i).p;
        end
        Seed_temp = [];
        for i = 1:length(NewSeedSet)
            CrowdDis = Crowding(NewSeedSet(i).p.decs);
            [~,R] = sort(CrowdDis);
            Seed_temp = [Seed_temp NewSeedSet(i).p(R(1))];
        end
        [SeedFrontNo,SeedMaxFNo] = NDSort(Seed_temp.objs,inf);
        while Sp_size < 2*N
            u = 0;
%             sigma = 0.04 * (SeedFrontNo(min_Sp)-1/SeedMaxFNo) + 0.01;
%             nor_obj = (Species(min_Sp).p.objs-min(Species(min_Sp).p.objs))./(max(Species(min_Sp).p.objs)-min(Species(min_Sp).p.objs));
%             sum_obj = sum(nor_obj,2);
%             [~,R] = sort(sum_obj);
%             ParentDec = Species(min_Sp).p(R(1)).dec;
            
            CrowdDis = Crowding(NewSeedSet(min_Sp).p.decs);
            [~,R] = sort(CrowdDis);
            ParentDec = NewSeedSet(min_Sp).p(R(1)).dec;
            [~,D]     = size(ParentDec);

%             a = 1/(10^randi([0,10])) .* (randi([1,9])+1);
%             b = ParentDec .* (a - ParentDec);
%             r = ParentDec ./ (a - ParentDec);
%             gaussian  = sigma * randn(1,D) + u;
%             OffspringDec = ParentDec + gaussian;
%             x=rand(1);%随机产生一个0到的数字
%             if x<0.5%如果小于0.5
%                 y=-1;%将-1付值给y
%             else
%                 y=1;%如果大于0.5，将1付值给y
%             end
            sigma = 0.01;
            gaussian  = sigma * randn(1,D)+u;
            OffspringDec = ParentDec + gaussian;
            if SeedFrontNo(i) > 1
%                 sigma = 0.01*SeedMaxFNo * (SeedFrontNo(i)/SeedMaxFNo);
                sigma = 0.01*SeedMaxFNo;
%                 sigma = 0.02 * (SeedFrontNo(i)-1);
%                   sigma = 0.05;
%                 z = rand(1);
%                 alpha = z*(1-exp(-0.9*Global.evaluation));
%                 if z > alpha
%                     sigma = 1-10^((-Global.D*Global.evaluated)/Global.evaluation);
%                 end
                gaussian  = sigma * randn(1,D) + u;
                OffspringDec = ParentDec + gaussian;
            end
            if SeedFrontNo(i) <= 1
                sigma = 0.01;
                gaussian  = sigma * randn(1,D) + u;
                OffspringDec = ParentDec + gaussian;
            end
            NewSeed = INDIVIDUAL(OffspringDec);
            Species(min_Sp).p = [Species(min_Sp).p, NewSeed];
            [FrontNo,MaxFNo] = NDSort([NewSeed.objs; NewSeedSet(min_Sp).p.objs],inf);
            if FrontNo(1) < FrontNo(2)
                NewSeedSet(min_Sp).p = NewSeed;
            elseif FrontNo(1) == FrontNo(2)
                NewSeedSet(min_Sp).p = [NewSeedSet(min_Sp).p, NewSeed];
            end
            length_Species = [];
            for i = 1:length(Species)
                length_Species = [length_Species length(Species(i).p)];
            end
            [~, max_Sp] = max(length_Species);
            [~, min_Sp] = min(length_Species);
            Sp_size = Sp_size + 1;
        end
        Population = [];
        for i = 1:length(Species)
            Population = [Population Species(i).p];
        end
        Seed = [];
%         for i = 1:length(Species)
%             [FrontNo,MaxFNo] = NDSort(Species(i).p.objs,inf);
%             First = find(FrontNo == 1);
%             Seed = [Seed Species(i).p(First)];
%         end
        for i = 1:length(NewSeedSet)
            Seed = [Seed NewSeedSet(i).p];
        end
        [FrontNo,MaxFNo] = NDSort(Seed.objs,inf);
        Seed_temp = [];
        for i = 1:length(Species)
            nor_obj = (Species(i).p.objs-min(Species(i).p.objs))./(max(Species(i).p.objs)-min(Species(i).p.objs));
            sum_obj = sum(nor_obj,2);
            [~,R] = sort(sum_obj);

%             CrowdDis = Crowding(Species(i).p.decs);
%             [~,R] = sort(CrowdDis);

            Seed_temp = [Seed_temp Species(i).p(R(1))];
        end
        [FrontNo,MaxFNo] = NDSort(Seed_temp.objs,inf);
        tabulate(FrontNo(:))
    end
     
    if mod == 2
        SeedSet = [];
        for i = 1:length(Species)
%             for j = 1:length(Species(i).p)
            nor_obj = (Species(i).p.objs-min(Species(i).p.objs))./(max(Species(i).p.objs)-min(Species(i).p.objs));
            sum_obj = sum(nor_obj,2);
            [~,R] = sort(sum_obj);
            SeedSet = [SeedSet Species(i).p(R(1))];
        end
        [FrontNo,MaxFNo] = NDSort(SeedSet.objs,inf);
%         [FrontNo,MaxFNo] = NDSort(center.objs,inf);
        tabulate(FrontNo(:))
        Seed = [];
        for i = 1:MaxFNo
            Front(i).p = [];
            Seed(i).p = [];
            for j = 1:length(Species)
                if FrontNo(j) == i
%                     for k = 1:length(Front(i-1))
%                     end
                    Front(i).p = [Front(i).p Species(j).p];
%                     [SpeciesFrontNo,~] = NDSort(Species(j).p.objs,inf);
%                     First = find(SpeciesFrontNo==1);
                    if length(Species(j).p) > 2
                        Seed(i).p = [Seed(i).p SeedSet(j)];
%                         if FrontNo(j) > 1
%                             [SpeciesFrontNo,~] = NDSort(Species(j).p.objs,inf);
%                              First = find(SpeciesFrontNo==1);
%                              Seed(i).p = [Seed(i).p Species(j).p(First)];
%                         end
                    end
%                     
                end
            end
        end
%         Seed = [];
%         for i = 1:length(Species)
% %             if length(Species(i).p) > 5
%                 [FrontNo,MaxFNo] = NDSort(Species(i).p.objs,inf);
%                 First = find(FrontNo==1);
%                 if length(Species(i).p(First)) > 5
%                     Seed = [Seed Species(i).p(First)];
%                 end
% %             end
% %             Seed = [Seed Species(i).p(First)];
%         end
       
%         %% 若此时物种规模小于N，则对选择最少物种中的目标空间拥挤距离最小值个体（种子）进行高斯分布的生成
%         while Sp_size < N
%     %         CrowdDis = Crowding(Population(Species(min_Sp).p).objs);
%     %         [~,R] = sort(CrowdDis,'descend');
%     %         ParentDec = Species(min_Sp).p(R(1)); %最后非支配层剔除位置索引
%             
%             
%             sigma = 0.1;
%     %         ParentDec = Species(i).p;
%             ParentDec = Seed(min_Sp).dec;
%             [~,D]     = size(ParentDec);
%             gaussian  = sigma * randn(1,D) + u;
%             OffspringDec = ParentDec + gaussian;
%             NewSeed = INDIVIDUAL(OffspringDec);
%     %         NewSeed = Global.Variation(Species(i).p,1,@Gaussian);
%             Population = [Population, NewSeed];
%             Species(min_Sp).p = [Species(min_Sp).p, NewSeed];
%             Seed = [Seed, NewSeed];
%             Sp_size = Sp_size + 1;
%             length_Species = [];
%             for i = 1:length(Species)
%                 length_Species = [length_Species length(Species(i).p)];
%             end
%             [m, index] = min(length_Species);
%             min_Sp = index;
%             min_Sp_size = m;
%         end
%         %% 若此时物种规模大于N，则对选择最大物种中的决策空间拥挤距离最小值个体进行剔除
%         while Sp_size > N
%             CrowdDis = Crowding(Species(max_Sp).p.decs);
%             [~,R] = sort(CrowdDis);
%     %         Population(Species(max_Sp).p(R(1))) = [];
%             Species(max_Sp).p(R(1)) = [];
%             Sp_size = Sp_size - 1;
%             length_Species = [];
%             for i = 1:length(Species)
%                 length_Species = [length_Species length(Species(i).p)];
%             end
%             [m, index] = max(length_Species);
%             max_Sp = index;
%             max_Sp_size = m;
%         end
        AllSpecies = [];
        for i = 1:length(Species)
            AllSpecies = [AllSpecies Species(i).p];
        end
    end
%     Seed = Population(Seed);
%     NewPopulation = Population(AllSpecies);
    
end