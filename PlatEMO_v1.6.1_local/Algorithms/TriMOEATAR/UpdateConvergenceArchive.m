function [AC,Rank,fS] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Global)
% Update the Convergence Archive

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

    Population = [Population,AC]; %将收敛性存储集与种群结合
    PopObj = Population.objs;
    PopDec = Population.decs;
    N = size(PopObj,1);   
    Rank   = inf(1,N);      % Rank of each solution
    nrank  = 1;             % Current rank
    
    %% Normalize
    PopObj = PopObj - repmat(Z,N,1);%对种群每个个体目标值减去最小值，进行标准化
    PopDec = (PopDec - repmat(Global.lower,N,1))./repmat(Global.upper-Global.lower,N,1);
    %对种群每个个体决策空间向量标准化，使其在[0,1]之间
    
    %% Convergence indicator
    fS = mean(PopObj'); %对每个个体计算其目标函数平均值 1*N矩阵
    % 每个解的收敛性指标
    %% Calculate distance between every two solutions in the IC decsion subspace
    d = pdist2(PopDec(:,Xic),PopDec(:,Xic),'chebychev');
    %收敛性解子空间每两个点的坐标数值差绝对值的最大值 切比雪夫距离
    %% Rank
    Choose = false(1,N);
    Q = true(1,N);
    Q1 = false(1,N);
    while sum(Choose) < NC %被选择的收敛集的规模没达到预设值
        if sum(Q) == 0 %若Q为空集 即最小收敛指标和其周围个体均被移除
           Q = Q1; %将未被直接选择但被划入nichies的个体可以重新选择
           Q1 = false(1,N);
           nrank = nrank+1;
        end
        % Choose x with min FS 
        temp1 = fS == min(fS(Q)); %在Q中找到收敛性指标最小的解的选择索引
        xmin = find(and(temp1,Q)); %所有最小解的位置索引
        xmin = xmin(1); %取第一个最小解
        Rank(xmin) = nrank; %其Rank值为当前nrank值
        Choose(xmin) = true; %置入Ac
        Q(xmin) = false; %从Q中删除
        % Delete solution near x_min
        temp3=d(xmin,:); %将xmin决策空间周围小于sigma半径的值 且在Q中
        temp2=temp3<sigma_niche;
        Delete = and(temp2,Q);
        Q(Delete) = false; %剔除Q
        Q1(Delete) = true;  %放入Q'
    end
    AC = Population(Choose);
    Rank = Rank(Choose);
    fS = fS(Choose);
end