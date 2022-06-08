function SCEAold(Global)
% <algorithm> <H-N>
% Fast and Elitist Multiobjective Genetic Algorithm: NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------
    % 200*N_ops  10000*N_ops
    Nns = Global.ParameterSet(Global.N/40);
    %% Generate random population
    [Z,Global.N] = UniformPoint(Global.N,Global.M);
    [ZD,Global.N] = UniformPoint(Global.N,Global.D);
    Population = Global.Initialization();
    Front_Zmin = [];
    Front_ZDmin = [];
    Zmin = min(Population.objs,[],1);
%     [~,FrontNo,CrowdDis] = EnvironmentalSelection(Population,Global.N);

    %% Optimization
    while Global.NotTermination(Population)
%         MatingPool = TournamentSelection(2,Global.N,FrontNo,-CrowdDis);
        if Global.evaluated < 0.5*Global.evaluation
            mod = 1;
        else
            mod = 2;
        end
        
        [Population,Species,center] = Species_division(Population,Nns,mod, Global);
        [Population, Seed, Species, Front] = Seed_detection(Population,Species,Global.N,mod, Global,center);
        if mod == 1
            Offspring = [];
            Zmin = min(Population.objs,[],1);
        end
        if mod == 2
            size_offspring = 0;
            for i = 1:length(Front)
                [FrontNo, ~] = NDSort(Front(i).p.objs,Front(i).p.cons,length(Front(i).p));
%                 MatingPool = randi(length(Front(i).p),1,floor(Global.N/length(Front)));
%                 MatingPool = TournamentSelection(2,length(Front(i).p),FrontNo);
                allFront = [];
%                 for j = 1:i
%                     allFront = [allFront Front(i).p];
%                 end
                allFront =  Front(i).p;
                MatingPool = randi(length(allFront),1,length(Front(i).p));
                Offspring  = Global.Variation(allFront(MatingPool));
                
%                     % MMODE生成新的子代
%                 MatingPool = randi(length(allFront),1,5*length(Front(i).p));
%                 Offspring  = Global.Variation(allFront(MatingPool),length(Front(i).p));

%                 MatingPool = randi(length(Front(i).p),1,length(Front(i).p));
%                 Offspring  = Global.Variation(Front(i).p(MatingPool));
%                 Front(i).p = [Front(i).p  Offspring(1:floor(Global.N/length(Front)))];
                Front(i).p = [Front(i).p  Offspring(1:length(Front(i).p))];
%                 size_offspring = size_offspring + floor(Global.N/length(Front));
                size_offspring = size_offspring + length(Front(i).p);
                Front_Zmin(i).p = min(Front(i).p.objs,[],1);
                Front_ZDmin(i).p = min(Front(i).p.decs,[],1);
            end
            if size_offspring < Global.N
                MatingPool = randi(length(Front(i).p),1,Global.N-size_offspring);
                Offspring  = Global.Variation(Front(i).p(MatingPool));
                Front(i).p = [Front(i).p  Offspring(1:Global.N-size_offspring)];
                Front_Zmin(i).p = min([Front_Zmin(i).p; Offspring(1:Global.N-size_offspring).objs],[],1);
                Front_ZDmin(i).p = min([Front_ZDmin(i).p; Offspring(1:Global.N-size_offspring).decs],[],1);
            end
        end
        Population = Seed_conservation([Population,Offspring],Seed,Species,Global.N,mod,Front,Global,Z,ZD,Zmin,Front_Zmin,Front_ZDmin,center);
    end
end