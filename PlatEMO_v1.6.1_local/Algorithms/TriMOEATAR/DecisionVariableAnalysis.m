function [Xic,Xre] = DecisionVariableAnalysis(Global,NCA,NIA)
% Decision variable analysis
% This code is modified from ControlVariableAnalysis.m and
% DividingDistanceVariables.m in MOEA-DVA

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
   
    Xic  = false(1,Global.D); %独立收敛性变量
    Xre  = false(1,Global.D); %
    
    %% Find convergence-related variable
    for i = 1 : Global.D
        x      = rand(1,Global.D).*(Global.upper-Global.lower) + Global.lower;
        S      = repmat(x,NCA,1); %NCA*D矩阵
        S(:,i) = ((1:NCA)'-1+rand(NCA,1))/NCA*(Global.upper(i)-Global.lower(i)) + Global.lower(i);
        [~,S]  = Global.problem('value',Global,S); %求初始化后改变第i个变量后所得目标值
        S = unique(S,'rows'); % delete the duplicate 删除重复项，按排序顺序返回S不重复的目标值
        [~,MaxFNo] = NDSort(S,inf);
        if MaxFNo == size(S,1) %非支配层数改变了
            Xic(i) = true; %则为独立收敛性相关
        else
            Xre(i) = true; %否则为相关性变量
        end
    end
    
    %% Interdependence 相关性
    % Generate the initial population        
    PopDec = rand(Global.N,Global.D);
    PopDec = PopDec.*repmat(Global.upper-Global.lower,Global.N,1) + repmat(Global.lower,Global.N,1); 
    %随机生成一个个体，并计算其目标值
    [~,PopObj]= Global.problem('value',Global,PopDec);
    % Interdependence analysis
    interaction = false(Global.D); %初始化关联逻辑矩阵
    interaction(logical(eye(Global.D))) = true; %单位矩阵转化为逻辑值，每个变量与其自身相关
    for i = 1 : Global.D-1
        for j = i+1 : Global.D
            for time2try = 1 : NIA
                % Detect whether the i-th and j-th decision variables are interacting
                x    = randi(Global.N);%1~N之间的随机整数
                a2   = rand*(Global.upper(i)-Global.lower(i)) + Global.lower(i); %随机取第i个变量
                b2   = rand*(Global.upper(j)-Global.lower(j)) + Global.lower(j); %随机取第j个变量
                Decs = repmat(PopDec(x,:),3,1); %随机取一个个体的决策空间向量，复制三份
                Decs(1,i) = a2;
                Decs(2,j) = b2;
                Decs(3,[i,j]) = [a2,b2]; %将其第i个、第j个和两个解向量分别替换作比较
                [~,F]= Global.problem('value',Global,Decs); %分别得到三组目标值
                delta1 = F(1,:) - PopObj(x,:); %j值不变仅改变第i个决策变量后的每个目标值的变化量
                delta2 = F(3,:) - F(2,:); %对另一个相同j值改变第i个决策变量
                interaction(i,j) = interaction(i,j) | any(delta1.*delta2<0); 
                %变化后的两组目标值中有任意目标值变化程度相反则为相关
                interaction(j,i) = interaction(i,j);                
            end
        end
    end
    
    %% Group based on Interdependence
    while sum(sum(interaction(Xic,Xre)))
        for i = find(Xic==1) %找到一个决策变量是独立的
            fprintf('i=%d\n',i);
            if sum(interaction(i,Xre)) %若它与某个变量相关
                Xic(i) = false; %则它并不独立
                Xre(i) = true; %放入相关联变量存储集中
            end
        end      
    end

end