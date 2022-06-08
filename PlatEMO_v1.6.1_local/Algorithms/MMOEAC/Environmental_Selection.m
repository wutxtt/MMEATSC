function Population=Environmental_Selection(Union, N, delta)
% The environmental selection
    
    %% Neighborhood-based Clustering Method
    range=max(Union.decs,[],1)-min(Union.decs,[],1);
    r=(range)*0.1; %取阈值r
    
    C=NCM(Union.decs,r); %决策空间的相邻解NCM聚类
    K=length(C); %聚类的个数
    local_A=[];
    % local_A: reserve the nondominated solutions in each local clusters
    for i=1:K
        cluster=C(i).p;
        if length(cluster)<=delta %预先设定好的在每个局部聚类中最少的解的个数
           continue;
        end
        [FrontNo1,~]=NDSort(Union(cluster).objs,length(cluster));
        local_A=[local_A cluster(FrontNo1==1)]; %将所有非支配解放入存储集中
    end
    
    %% P: the number of individuals is not less than N
    [FrontNo,~]=NDSort(Union.objs,length(Union));
    P = find(FrontNo==1);   %对整个种群非支配排序
    temp=setdiff(P,local_A);
    P=[local_A temp];
    i=1;
    while  length(P)<N
        i=i+1;
        temp1=find(FrontNo==i);
        temp2=setdiff(temp1, local_A);
        P=[P temp2]; 
    end %若种群数量不足N则添加非重复的正常的非支配解

    if length(P)>N 
        newpop=Union(P);
        z=min(newpop.objs,[],1);
        Z=max(newpop.objs,[],1);
        % normalization in the decision space 对目标空间的标准化
        ZZZ=(newpop.objs-z)./repmat(Z-z,length(newpop),1);
        % hierarchical clustering method
        H=clusterdata(ZZZ,'maxclust',N,'distance','euclidean','linkage','ward');
        % ward's linkage HCM目标空间的聚类 若聚类簇的数目大于种群规模N则继续聚类
        count=length(newpop)-N;
        %% delete redundant individual in turn
        for i=1:count %删除种群中大于N数目的个体
            CrowdDis=Crowding(newpop.decs);%计算每个解决策空间的HAD
            num=hist(H,1:N);
            % find the most crowding clusters in the objective space 
            % and then add them into R
            I=find(num==max(num)); %即解最多的聚类
            R=[];
            for j=1:length(I)
                R=[R;find(H==I(j))];
            end
            % delete the solution with minimal decision spatial crowding distance 
            T=find(min(CrowdDis(R))==CrowdDis(R)); %有最小HAD的位置索引
            s=randperm(length(T)); %随机选择有最小HAD的解位置索引
            x=R(T(s(1))); %找到选择解的对应索引 在种群和簇中删除
            newpop(x)=[];
            H(x)= [];
        end
        Population=newpop;
    else
        Population=Union(P);
    end
end