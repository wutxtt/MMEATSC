function Population = EnvironmentalSelection(Population, Global, NS)
% The environmental selection of NSGA-II

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    
    for i = 1:length(Population)
        NP = length(Population); %种群规模
        decc = Population.decs;
        d_dec = pdist2(decc,decc,'euclidean');  %计算决策空间每个点之间的距离
        NB = min(NS,max(NP-1,0)); %取有最小距离个体的个数
        temp = sort(d_dec); %对于Nns个个体之间的距离进行按列升序排序
        r_dec = sum(temp(i,1:NB))/NB; %计算NB*NP个距离的平均距离
        neighbor = [];
        for j = length(d_dec(i,:))
            if d_dec(i,j) > 0 && d_dec(i,j) < r_dec
                neighbor = [neighbor j];
            end
        end
        NEI = [Population(i),Population(neighbor)];
        PCD = 0;
        for j = 1:length(NEI)
            p = Population(i).decs;
            n = NEI(j).decs;
            PCD = PCD + norm(p-n);
        end
        PCD = PCD / length(NEI);
        [FrontNo,MaxFNo] = NDSort(NEI.objs,NEI.cons,length(NEI));
        First = FrontNo == 1;
        PFno = FrontNo(1);
        F = NEI(First);
        Offspring  = Global.Variation([F,Population(i)],1);
        
        OCD = 0;
        for j = 1:length(NEI)
            OCD = OCD + norm(Offspring.decs-NEI(j).decs);
        end
        OCD = OCD / length(NEI);
        [FrontNo,MaxFNo] = NDSort([Offspring.objs;NEI.objs],[Offspring.cons;NEI.cons],length(NEI)+1);
        First = FrontNo == 1;
        OFno = FrontNo(1);
        
        if OFno < PFno 
            Population(i) = Offspring;
        end
        if OFno == PFno &&  OCD > PCD
            Population(i) = Offspring;
        end
    end
end