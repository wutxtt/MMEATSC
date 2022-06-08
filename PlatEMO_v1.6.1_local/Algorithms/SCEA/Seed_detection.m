function [Population, Seed, Species, Front] = Seed_detection(Population,Species,N,mod, Global,center, r_dec)
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
            if 
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
%         Seed_t = [];
%         for i = 1:length(Seed)
%             Seed_t = [Seed_t Seed(i).p];
%         end
%         Seed = Seed_t;
        [SeedFrontNo,SeedMaxFNo] = NDSort(Seed_temp.objs,inf);
        while Sp_size < 2*N
            u = 0;
            CrowdDis = Crowding(NewSeedSet(min_Sp).p.decs);
            [~,R] = sort(CrowdDis);
            ParentDec = NewSeedSet(min_Sp).p(R(1)).dec;
            [~,D]     = size(ParentDec);
%             sigma = 0.001;
%             gaussian  = sigma * randn(1,D)+u;
%             OffspringDec = ParentDec + gaussian;
            if SeedFrontNo(i) > 3
                sigma = 0.01*SeedFrontNo(i);
                gaussian  = sigma * randn(1,D) + u;
                OffspringDec = ParentDec + gaussian;
            end
            if SeedFrontNo(i) <= 3
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
        Seed = [];
        for i = 1:length(Species)
            Population = [Population Species(i).p];
%             [FrontNo,MaxFNo] = NDSort(Species(i).p.objs,Species(i).p.cons,length(Species(i).p));
%             First = find(FrontNo==1);
%             Seed = [Seed Species(i).p(First)];
        end
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
%         Seed = Seed_temp;
    end
     
    if mod == 2
        
        all_seed_set = [];
        for i = 1:length(Species)
            Population = [Population Species(i).p];
            [FrontNo,MaxFNo] = NDSort(Species(i).p.objs,Species(i).p.cons,length(Species(i).p));
            First = find(FrontNo==1);
            all_seed_set(i).p = Species(i).p(First);
        end
        
        SeedSet = [];
        for i = 1:length(Species)
            nor_obj = (Species(i).p.objs-min(Species(i).p.objs))./(max(Species(i).p.objs)-min(Species(i).p.objs));
            sum_obj = sum(nor_obj,2);
            [~,R] = sort(sum_obj);
            SeedSet = [SeedSet Species(i).p(R(1))];
        end
        [FrontNo,MaxFNo] = NDSort(SeedSet.objs,inf);
%         [FrontNo,MaxFNo] = NDSort(center.objs,inf);
        tabulate(FrontNo(:))
        Seed = [];
        f = 0; %f为当前分组下标
        r = 2*r_dec;
        b = 5;
        LastFront = [];
        LastFrontSeed = [];
        for i = 1:MaxFNo %i为帕累托等级
%             if f == 0
%                 f = f + 1;
%                 FrontSeed(f).p = [];
%             end
%             if i > 4 && i > MaxFNo/2
             if i > 3
% %                 for j = 1:length(Species) %j为物种下标
% %                     if FrontNo(j) == i
% %                         Front(f).p = [Front(f).p Species(j).p];
% %                         FrontSeed(f).p = [FrontSeed(f).p SeedSet(j)];
% %                         if length(Species(j).p) > 2
% %                             Seed(f).p = [Seed(f).p SeedSet(j)];
% %                         end
% %                     end
% %                 end
% %                 temp_f = f; % temp_f为当前帕累托等级的第一个分组下标

%                 for j = 1:length(Species) %j为物种下标
%                     if FrontNo(j) == i
%                         LastFront = [LastFront Species(j).p];
%                     end
%                 end

                if f == 0
                    f = f + 1;
                    FrontSeed(f).p = [];
                    temp_f = f;
                end
                for j = 1:length(Species) %j为物种下标
                    if FrontNo(j) == i
                        if isempty(FrontSeed(f).p)
                            Front(f).p = Species(j).p;
                            FrontSeed(f).p = SeedSet(j);
                            Seed(f).p = [];
%                             if length(Species(j).p) > b
%                                 Seed(f).p = SeedSet(j);
% %                                 Seed(f).p = all_seed_set(j).p;
%                             end
                        else
                            flag = 0;
                            for k = temp_f:f
                                tempseed = FrontSeed(k).p;
                                for l = 1:length(tempseed) %l为当前分组的种子下标
                                    Distance = pdist([tempseed(l).decs;SeedSet(j).decs],'euclidean');
                                    if Distance < r
                                        Front(k).p = [Front(k).p Species(j).p];
                                        FrontSeed(k).p = [FrontSeed(k).p SeedSet(j)];
%                                         if length(Species(j).p) > b
%                                             Seed(k).p = [Seed(k).p SeedSet(j)];
% %                                             Seed(f).p = [Seed(k).p all_seed_set(j).p];
%                                         end
                                        flag = 1;
                                        break
                                    end
                                end
                                if flag == 1
                                    break
                                end
                            end
                            if flag == 0 
                                f = f + 1;
                                Front(f).p = Species(j).p;
                                FrontSeed(f).p = SeedSet(j);
                                Seed(f).p = [];
%                                 if length(Species(j).p) > b
%                                     Seed(f).p = SeedSet(j);
% %                                     Seed(f).p = all_seed_set(j).p;
%                                 end
                            end 
                        end     
                    end
                end
            else
%                 if f ~= 1
                    f = f + 1;
%                 end
                FrontSeed(f).p = [];
                temp_f = f; % temp_f为当前帕累托等级的第一个分组下标
                for j = 1:length(Species) %j为物种下标
                    if FrontNo(j) == i
                        if isempty(FrontSeed(f).p)
                            Front(f).p = Species(j).p;
                            FrontSeed(f).p = SeedSet(j);
                            Seed(f).p = [];
                            if length(Species(j).p) > b
                               Seed(f).p = SeedSet(j);
%                                 Seed(f).p = all_seed_set(j).p;
                            end
                        else
                            flag = 0;
                            for k = temp_f:f
                                tempseed = FrontSeed(k).p;
                                for l = 1:length(tempseed) %l为当前分组的种子下标
                                    Distance = pdist([tempseed(l).decs;SeedSet(j).decs],'euclidean');
                                    if Distance < r
                                        Front(k).p = [Front(k).p Species(j).p];
                                        FrontSeed(k).p = [FrontSeed(k).p SeedSet(j)];
                                        if length(Species(j).p) > b
                                            Seed(k).p = [Seed(k).p SeedSet(j)];
%                                             Seed(f).p = [Seed(k).p all_seed_set(j).p];
                                        end
                                        flag = 1;
                                        break
                                    end
                                end
                                if flag == 1
                                    break
                                end
                            end
                            if flag == 0 
                                f = f + 1;
                                Front(f).p = Species(j).p;
                                FrontSeed(f).p = SeedSet(j);
                                Seed(f).p = [];
                                if length(Species(j).p) > b
                                    Seed(f).p = SeedSet(j);
%                                     Seed(f).p = all_seed_set(j).p;
                                end
                            end 
                        end     
                    end
                end
            end
        end
        
%         for i = 1:MaxFNo
%             Front(i).p = [];
%             Seed(i).p = [];
%             for j = 1:length(Species)
%                 if FrontNo(j) == i
%                     Front(i).p = [Front(i).p Species(j).p];
%                     if length(Species(j).p) > 2
%                         Seed(i).p = [Seed(i).p SeedSet(j)];
%                     end          
%                 end
%             end
%         end
%         if ~isempty(LastFront)
%                Front(f+1).p = LastFront;
%                FrontSeed(f+1).p = [];
%         end
%        
        AllSpecies = [];
        for i = 1:length(Species)
            AllSpecies = [AllSpecies Species(i).p];
        end
    end
%     Seed = Population(Seed);
%     NewPopulation = Population(AllSpecies);
    
end