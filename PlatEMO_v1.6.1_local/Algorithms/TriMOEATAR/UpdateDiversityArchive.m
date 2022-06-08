function [AD,Rank] = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Global)
% Update the Diversity Archive

%--------------------------------------------------------------------------
% Copyright 2017-2018 Yiping Liu
% This is the code of TriMOEA-TA&R proposed in "Yiping Liu, Gary G. Yen, 
% and Dunwei Gong, A Multi-Modal Multi-Objective Evolutionary Algorithm 
% Using Two-Archive and Recombination Strategies, IEEE Transactions on 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    %% Non-dominated sorting
    Population = [Population,AD];
    [FrontNo,MaxFNo] = NDSort(Population.objs,ND); %ND为多样性存储集的规模

    %% Select the solutions in the last front
    if length(find(FrontNo<=MaxFNo))== ND
       Next = FrontNo<=MaxFNo;
    else
       Next = FrontNo < MaxFNo ;
       NQ     = ND-sum(Next);
       Last   = find(FrontNo==MaxFNo);
       Choose = LastSelection(Population(Last).objs,Population(Last).decs,NQ,R,Z,Xre,sigma_niche,Global); %Algorithm 4, lines 5-30
       Next(Last(Choose)) = true;
    end
    
    %% Population for next generation
    AD = Population(Next);
    Rank = FrontNo(Next);
end

function Choose = LastSelection(PopObj,PopDec,NQ,R,Z,Xre,sigma_niche,Global)
% Select solutions with good diversity in the last front    
    N      = size(PopObj,1);
    NR     = size(R,1);
    
    %% Normalize
    PopObj = PopObj - repmat(Z,N,1);
    PopDec = (PopDec - repmat(Global.lower,N,1))./repmat(Global.upper-Global.lower,N,1);
  
    %% Calculate distance between every solution and reference point in the objective space
    theta = pdist2(R,PopObj,'cosine');
    %计算参考点与标准化目标值之间的角度距离 夹角越小余弦值越小
    %% Calculate distance between every two solutions in the reminder decsion subspace
    d = pdist2(PopDec(:,Xre),PopDec(:,Xre),'chebychev');   
    %计算其余解在剩余子空间的切比雪夫距离
    %% Cluster
    [thmin,label] = min(theta,[],1);   %返回1*N 每个个体相关的参考点 及参考点对应的索引
    C1 = false(NR,N); %每个参考点的相关个体
    C2 = false(NR,N);
    for j = 1:NR
        member = find(label==j); %该参考点相关个体的位置索引
        member1 = label==j; %该参考点相关个体的选择索引
        [~,temp] = sort(thmin(member)); %对所有相关该参考点的个体进行余弦值升序排序
        for i = temp
           if any(d(member(i),and(C1(j,:),member1))<sigma_niche) 
               C2(j,member(i))=true;
           else
               C1(j,member(i))=true;
           end
        end        
    end    
    
    %% Make selected solution == NQ
    while sum(sum(C1))>NQ
        cmax = max(sum(C1,2));
        jmax = sum(C1,2)== cmax;
        temp1 = find(sum(C1(jmax,:),1)>0);        
        [~,xmax] = max(thmin(temp1));
        xmax=temp1(xmax);
        C1(label(xmax),xmax) = false;
    end   
    while sum(sum(C1))<NQ
        c2 = sum(C2,2)>0;
        cmin = min(sum(C1(c2,:),2));
        jmin = and(sum(C1,2)==cmin,c2);
        temp1 = find(sum(C2(jmin,:),1)>0);
        [~,xmin] = min(thmin(temp1));
        xmin=temp1(xmin);
        C2(label(xmin),xmin) = false;
        C1(label(xmin),xmin) = true;
    end    
    Choose = sum(C1)>0;  
end