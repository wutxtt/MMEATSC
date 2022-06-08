function TriMOEATAR(Global)
% <algorithm> <MMMO>
% A Multi-Modal Multi-Objective Evolutionary Algorithm Using Two-Archive 
% and Recombination Strategies
% p_con --- 0.5 --- Probability of selecting parents from the convergence archive
% sigma_niche --- 0.1 --- Niche radius in the decision space
% eps_peak --- 0.01 --- Accuracy level to detect peaks
% NR --- 100 --- The number of refernece points
% NCA --- 20 --- The number of sampling solutions in control variable analysis
% NIA --- 6 --- The maximum number of tries required to judge the interaction

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

    %% Parameter setting
    [p_con,sigma_niche,eps_peak,NR,NCA,NIA] = Global.ParameterSet(0.5,0.1,0.01,100,20,6);  
    %选择收敛性解的占比 变量空间的小生境半径 融合策略的准确度 参考点的个数 判断独立收敛性的样本选择数 判断相关性的变量改变次数
    %% Decision Variable Analysis and Generate random population
    [Xic,Xre] = DecisionVariableAnalysis(Global,NCA,NIA);
    
    %% Generate the reference points    
    [R,~] = UniformPoint(NR,Global.M);         
    %在M维目标空间中，返回NR个均匀分布的参考点
    %% Initialization
    Population = Global.Initialization();
    NC = Global.N;
    ND = Global.N;
    AC = [];
    AD = [];
    Pc = floor(p_con.*Global.N); %选取的收敛性相关个体，向着负无穷舍入，将每个元素舍入到小于或等于该元素的最接近整数
    Pd = Global.N - Pc;      
    
    %% Update Archives
    Z  = min(Population.objs,[],1);  %对每列取最小值返回行，取所有目标值的最小值
    [AC,RankC,~] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Global);
    [AD,RankD] = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Global);

    %% Optimization
    while Global.NotTermination(Population)
        %% Generate new Population based on two Archives
        MatingPoolC = TournamentSelection(2,Pc,RankC);
        MatingPoolD = TournamentSelection(2,Pd,RankD);
        Parents     = [AC(MatingPoolC),AD(MatingPoolD)];
        Parents     = Parents(randperm(Global.N));
        Population  = Global.Variation(Parents);        
        %% Update Archives
        Z             = min([Z;Population.objs],[],1);
        [AC,RankC,fS] = UpdateConvergenceArchive(AC,Population,NC,Z,Xic,sigma_niche,Global);        
        [AD,RankD]    = UpdateDiversityArchive(AD,Population,ND,R,Z,Xre,sigma_niche,Global);                
        %% Recombination
        if Global.gen == Global.maxgen
            if sum(Xic) == 0
               Population = AD;
            else
               FS = Recombination(AC,AD,Xic,Xre,eps_peak,fS,RankC);
               Population = INDIVIDUAL(FS);
            end
        end        
    end    
end