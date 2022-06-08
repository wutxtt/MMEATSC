function Population = Clearing(Population,Offspring,Global)
    Population = [Population, Offspring];
    Ns = Global.N/(2*(Global.D-1)) + ((Global.N-Global.N/(2*(Global.D-1))) * sin((Global.evaluated/Global.evaluation)*pi));
    num = length(Population);

    dist = pdist2(Population.decs,Population.decs);
    dist = sort(dist);
    d = sort(dist);
    dd = d(:,1:Ns);
    dist = sum(sum(dd));

    r_dec = dist/(Ns*num);
%%
    Archive = [];
%     S = (1:num);  % S为1*n的矩阵
    d_dec = pdist2(Population.decs,Population.decs,'euclidean');  %计算决策空间每个点之间的距离
    for x = 1:num
        % C为包含全部物种的集合
        center = Population(x);
        C(x).p = [];  % 初始化一个物种中个体的list
        Distance = d_dec(x,:);  % 取这个中心个体与全部其他个体在决策空间的距离
        Nbs = find(Distance < r_dec);  % 得到所有在决策空间距离小于r_dec的个体
        C(x).p = [C(x).p Nbs];  % 将所得到的个体全部放入初始化的list中
        
        C(x).p = setdiff(C(x).p,Archive);
%         C(x).p = [x C(x).p];
        S = [center Population(C(x).p)];
        [FrontNo, ~] = NDSort(S.objs,length(C(x).p));
        if FrontNo(1) > 1
            Archive = [Archive x];
        end
    end
    
    %% 对A正则化
    A = Population(Archive);
    z=min(A.objs,[],1);
    Z=max(A.objs,[],1);
    
    nf_obj=(A.objs-z)./repmat(Z-z,length(A),1);
    
    %% 计算A中任意解的收敛性指标
    C = [];
    for i = 1:length(A)
        s = 0;
        for l = 1:length(A(i).objs)
            o = A(i).objs;
            s = s + o(l)^2;
        end
        s = sqrt(s);
        C = [C s];
    end
    P = Population;
    P(Archive) = [];
    %% 找到A中最大和最小收敛指标值的个体c1,c2
    while (length(A) + length(P)) > Global.N
        if length(A) == 0
            break
        end
        [~,c1] = min(C);
        A(c1) = [];
        C(c1) = [];
        [~,c2] = max(C);
        A(c2) = [];       
        C(c2) = [];
    end
    
    Population = [A P];

end