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
%     S = (1:num);  % SΪ1*n�ľ���
    d_dec = pdist2(Population.decs,Population.decs,'euclidean');  %������߿ռ�ÿ����֮��ľ���
    for x = 1:num
        % CΪ����ȫ�����ֵļ���
        center = Population(x);
        C(x).p = [];  % ��ʼ��һ�������и����list
        Distance = d_dec(x,:);  % ȡ������ĸ�����ȫ�����������ھ��߿ռ�ľ���
        Nbs = find(Distance < r_dec);  % �õ������ھ��߿ռ����С��r_dec�ĸ���
        C(x).p = [C(x).p Nbs];  % �����õ��ĸ���ȫ�������ʼ����list��
        
        C(x).p = setdiff(C(x).p,Archive);
%         C(x).p = [x C(x).p];
        S = [center Population(C(x).p)];
        [FrontNo, ~] = NDSort(S.objs,length(C(x).p));
        if FrontNo(1) > 1
            Archive = [Archive x];
        end
    end
    
    %% ��A����
    A = Population(Archive);
    z=min(A.objs,[],1);
    Z=max(A.objs,[],1);
    
    nf_obj=(A.objs-z)./repmat(Z-z,length(A),1);
    
    %% ����A��������������ָ��
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
    %% �ҵ�A��������С����ָ��ֵ�ĸ���c1,c2
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