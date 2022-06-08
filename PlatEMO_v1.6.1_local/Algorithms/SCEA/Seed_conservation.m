function [Population] = Seed_conservation(Population,Seed_population,Species,N,mod,Front,Global,Z,ZD,Zmin,ZDmin,Front_Zmin,Front_ZDmin,center)
% The environmental selection of NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
%     Seed_population = Population(Seed);
%     %% normalization in the decision space 对决策空间的标准化
%     z=min(Population.objs,[],1);
%     Z=max(Population.objs,[],1);
%     objj=(Population.objs-z)./repmat(Z-z,length(Population),1);
%     %% normalization in the decision space 对决策空间的标准化
%     z=min(Population.decs,[],1);
%     Z=max(Population.decs,[],1);
%     decc=(Population.decs-z)./repmat(Z-z,length(Population),1);
    
    %% FEs < 0.5*gen
    if mod == 1
%         Sp_size = size(Population,2);
%         CrowdDis = Crowding(Population.decs);
%         [~,R] = sort(CrowdDis);
%         Population(R(1:N+length(Seed_population))) = [];
        % Non-dominated sorting
        pobjs = Population.objs;
%         [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N-length(Seed_population));
        [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,N);
        Next = FrontNo < MaxFNo;

        % Calculate the crowding distance of each solution
%         CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        

        % Select the solutions in the last front based on their crowding distances
%         Last     = find(FrontNo==MaxFNo);
%         CrowdDis = Crowding(Population(Last).objs);
% %         [~,Rank] = sort(CrowdDis(Last),'descend');
%         [~,Rank] = sort(CrowdDis,'descend');
%         Next(Last(Rank(1:N-length(Seed_population)-sum(Next)))) = true;
% %         Next(Last(Rank(1:N-sum(Next)))) = true;

        % 基于NSGAIII分解参考点方法的选择
%         if length(Seed_population) > N
%             Population = 
        Last   = find(FrontNo==MaxFNo);
%         Choose = LastSelection(Population(Next).objs,Population(Last).objs,N-length(Seed_population)-sum(Next),Z,Zmin);
        Choose = LastSelection(Population(Next).decs,Population(Last).decs,N-length(Seed_population)-sum(Next),ZD,ZDmin);
%         Choose = LastSelection(Population(Next).decs,Population(Last).decs,N-sum(Next),ZD,ZDmin);
        Next(Last(Choose)) = true;
        
        % Population for next generation
%         pp = Population.decs;
%         clf
%         plot(pp(:,1),pp(:,2),'o');
        
        Population = Population(Next);
        
%         pp = Population.decs;
%         if length(pp) > 0
%             plot(pp(:,1),pp(:,2),'o');
%         else
%             clf
%         end
%         ss = Seed_population.decs;
%         hold on
%         plot(ss(:,1),ss(:,2),'ro');
        
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
%         Seed = [];
%         for i = 1:MaxFNo
%             Front(i).p = [];
%             Seed(i).p = [];
%             for j = 1:length(Species)
%                 if FrontNo(j) == i
%                     Front(i).p = [Front(i).p Species(j).p];
%                     if length(Species(j).p) > 2
%                         Seed(i).p = [Seed(i).p SeedSet(j)];
%                     end
% %                     
%                 end
%             end
%         end
        
%         for i = 1:length(Front)
%             ff = Front(i).p.decs;
%             hold on
%             plot(ff(:,1),ff(:,2),'*');
%         end

%         nor_obj = (Population.objs-min(Population.objs))./(max(Population.objs)-min(Population.objs));
%         sum_obj = sum(nor_obj,2);
%         [~,R] = sort(sum_obj);
%         Population(R(1:min(N,length(Seed_population)))) = [];
        
        Population = [Population, Seed_population];
    end

    %% FEs > 0.5*gen
    if mod == 2
        for i = 1:length(Front)
            [FrontNo,MaxFNo] = NDSort(Front(i).p.objs,Front(i).p.cons,length(Front(i).p)/2);
            Next = FrontNo < MaxFNo;
            
%             Select the solutions in the last front based on their crowding distances
% %             CrowdDis = CrowdingDistance(Front(i).p.objs,FrontNo);
%             Last     = find(FrontNo==MaxFNo);
% %             CrowdDis = Crowding(Front(i).p(Last).objs);
%             CrowdDis = Crowding(Front(i).p(Last).decs);
% %             [~,Rank] = sort(CrowdDis(Last),'descend');
%             [~,Rank] = sort(CrowdDis,'descend');
% %             Next(Last(Rank(1:length(Front(i).p)/2-length(SeedSet(i).p)-sum(Next)))) = true;
%             Next(Last(Rank(1:length(Front(i).p)/2-sum(Next)))) = true;
            
            
%             基于NSGAIII分解参考点方法的选择
            Last   = find(FrontNo==MaxFNo);
%             Choose = LastSelection(Front(i).p(Next).objs,Front(i).p(Last).objs,length(Front(i).p)/2-sum(Next),Z,Front_Zmin(i).p);
            Choose = LastSelection(Front(i).p(Next).decs,Front(i).p(Last).decs,length(Front(i).p)/2-sum(Next),ZD,Front_ZDmin(i).p);
%             Choose = LastSelection(Front(i).p(Next).decs,Front(i).p(Last).decs,length(Front(i).p)/2-sum(Next)-length(Seed_population(i).p),ZD,Front_ZDmin(i).p);
            Next(Last(Choose)) = true;
            
            Front(i).p = Front(i).p(Next);
            
%             Choose = LastSelection([],Front(i).p.decs,length(Front(i).p)/2,ZD,Front_ZDmin(i).p);
%             Front(i).p = Front(i).p(Choose);
            
        end
        Population = [];
        for i = 1:length(Front)
            Population = [Population Front(i).p];
        end
%         clf
%         pp = Population.decs;
%         plot(pp(:,1),pp(:,2),'o');
%         hold on
%         for i = 1:length(Seed_population)
%             ss = Seed_population(i).p.decs; 
%             plot(ss(:,1),ss(:,2),'ro');
%             hold on
%         end
        Seed = [];
        for i = 1:length(Seed_population)
            Seed = [Seed Seed_population(i).p];
        end
        if Global.evaluated < Global.evaluation
%             nor_obj = (Population.objs-min(Population.objs))./(max(Population.objs)-min(Population.objs));
%             sum_obj = sum(nor_obj,2);
%             [~,R] = sort(sum_obj);
%             Population(R(1:length(Seed_population))) = [];
%             CrowdDis = Crowding(Population.decs);
%             
%             [~,R] = sort(-CrowdDis);
%             Population(R(1:min(N,length(Seed)))) = [];
%             
%             Population = [Population, Seed_population];
            Population = [Population Seed];
        end
        
    end
    if mod == 3
        Population = [];
        Seed = [];
%         for i = 1:length(Seed_population)
%             Seed = [Seed Seed_population(i).p];
%         end
        for i = 1:length(Front)
            Population = [Population Front(i).p];
        end
        [FrontNo,MaxFNo] = NDSort(Population.objs,Population.cons,2*N-length(Seed_population));
        Next = FrontNo < MaxFNo;
        CrowdDis = CrowdingDistance(Population.objs,FrontNo);
        Last     = find(FrontNo==MaxFNo);
        [~,Rank] = sort(CrowdDis(Last),'descend');
        Next(Last(Rank(1:2*N-length(Seed_population)-sum(Next)))) = true;
        Population = Population(Next);
        Population = [Population, Seed_population];
    end
   
end


function Choose = LastSelection(PopObj1,PopObj2,K,Z,Zmin)
% Select part of the solutions in the last front

    PopObj = [PopObj1;PopObj2] - repmat(Zmin,size(PopObj1,1)+size(PopObj2,1),1);
    [N,M]  = size(PopObj);
    N1     = size(PopObj1,1);
    N2     = size(PopObj2,1);
    NZ     = size(Z,1);

    %% Normalization
    % Detect the extreme points
    Extreme = zeros(1,M);
    w       = zeros(M)+1e-6+eye(M);
    for i = 1 : M
        [~,Extreme(i)] = min(max(PopObj./repmat(w(i,:),N,1),[],2));
    end
    % Calculate the intercepts of the hyperplane constructed by the extreme
    % points and the axes
    Hyperplane = PopObj(Extreme,:)\ones(M,1);
    a = 1./Hyperplane;
    if any(isnan(a))
        a = max(PopObj,[],1)';
    end
    % Normalization
    PopObj = PopObj./repmat(a',N,1);
    
    %% Associate each solution with one reference point
    % Calculate the distance of each solution to each reference vector
    Cosine   = 1 - pdist2(PopObj,Z,'cosine');
    Distance = repmat(sqrt(sum(PopObj.^2,2)),1,NZ).*sqrt(1-Cosine.^2);
    % Associate each solution with its nearest reference point
    [d,pi] = min(Distance',[],1);

    %% Calculate the number of associated solutions except for the last front of each reference point
    rho = hist(pi(1:N1),1:NZ);
    
    %% Environmental selection
    Choose  = false(1,N2);
    Zchoose = true(1,NZ);
    % Select K solutions one by one
    while sum(Choose) < K
        % Select the least crowded reference point
        Temp = find(Zchoose);
        Jmin = find(rho(Temp)==min(rho(Temp)));
        j    = Temp(Jmin(randi(length(Jmin))));
        I    = find(Choose==0 & pi(N1+1:end)==j);
        % Then select one solution associated with this reference point
        if ~isempty(I)
            if rho(j) == 0
                [~,s] = min(d(N1+I));
            else
                s = randi(length(I));
            end
            Choose(I(s)) = true;
            rho(j) = rho(j) + 1;
        else
            Zchoose(j) = false;
        end
    end
end