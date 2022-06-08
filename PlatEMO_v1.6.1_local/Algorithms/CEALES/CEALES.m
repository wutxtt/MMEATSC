function CEALES(Global)
% <algorithm> <A>
% Hierarchy Ranking Based Multimodal Multi-objective Evolutionary Algorithm
% eps --- 0.3 --- parameter for quality of the local PF
% p --- 0.5 --- parameter for probability

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
[k] = Global.ParameterSet(5);

%% Generate random population
Population = Global.Initialization();


%% Optimization
while Global.NotTermination(Population)
    [FrontNo, MaxFNo] = NDSort(Population.objs,Global.N);
    for i = 1:MaxFNo
        temp = FrontNo == i;
        Subpop(i).p = Population(temp);
    end
    Offspring = [];
    for i = MaxFNo:-1:1
%         temp = [Subpop(i).p,Subpop(i-1).p];
        if i ~= 1
            temp = [Subpop(i).p,Subpop(i-1).p];
        else
            temp = [Subpop(i).p];
        end
        if Global.evaluated <= Global.evaluation * 0.5
            MatingPool = randi(length(Subpop(i).p),1,length(temp));
%             MatingPool = TournamentSelection(2,round(Global.N),-CrowdDis2);
            O = Global.Variation([temp(MatingPool)],length(Subpop(i).p));
        else
            [FrontNo, ~] = NDSort(temp.objs,Global.N);
            MatingPool = TournamentSelection(2,round(Global.N),FrontNo);
            O  = Global.Variation([temp(MatingPool)],length(Subpop(i).p));
        end
        Offspring = [Offspring, O];
    end
    Population = Clearing(Population, Offspring, Global);
    Population = EnvironmentalSelection(Population,Global.N,k);
end

end
