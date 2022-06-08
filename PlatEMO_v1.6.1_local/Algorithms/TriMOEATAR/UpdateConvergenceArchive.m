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

    Population = [Population,AC]; %�������Դ洢������Ⱥ���
    PopObj = Population.objs;
    PopDec = Population.decs;
    N = size(PopObj,1);   
    Rank   = inf(1,N);      % Rank of each solution
    nrank  = 1;             % Current rank
    
    %% Normalize
    PopObj = PopObj - repmat(Z,N,1);%����Ⱥÿ������Ŀ��ֵ��ȥ��Сֵ�����б�׼��
    PopDec = (PopDec - repmat(Global.lower,N,1))./repmat(Global.upper-Global.lower,N,1);
    %����Ⱥÿ��������߿ռ�������׼����ʹ����[0,1]֮��
    
    %% Convergence indicator
    fS = mean(PopObj'); %��ÿ�����������Ŀ�꺯��ƽ��ֵ 1*N����
    % ÿ�����������ָ��
    %% Calculate distance between every two solutions in the IC decsion subspace
    d = pdist2(PopDec(:,Xic),PopDec(:,Xic),'chebychev');
    %�����Խ��ӿռ�ÿ�������������ֵ�����ֵ�����ֵ �б�ѩ�����
    %% Rank
    Choose = false(1,N);
    Q = true(1,N);
    Q1 = false(1,N);
    while sum(Choose) < NC %��ѡ����������Ĺ�ģû�ﵽԤ��ֵ
        if sum(Q) == 0 %��QΪ�ռ� ����С����ָ�������Χ��������Ƴ�
           Q = Q1; %��δ��ֱ��ѡ�񵫱�����nichies�ĸ����������ѡ��
           Q1 = false(1,N);
           nrank = nrank+1;
        end
        % Choose x with min FS 
        temp1 = fS == min(fS(Q)); %��Q���ҵ�������ָ����С�Ľ��ѡ������
        xmin = find(and(temp1,Q)); %������С���λ������
        xmin = xmin(1); %ȡ��һ����С��
        Rank(xmin) = nrank; %��RankֵΪ��ǰnrankֵ
        Choose(xmin) = true; %����Ac
        Q(xmin) = false; %��Q��ɾ��
        % Delete solution near x_min
        temp3=d(xmin,:); %��xmin���߿ռ���ΧС��sigma�뾶��ֵ ����Q��
        temp2=temp3<sigma_niche;
        Delete = and(temp2,Q);
        Q(Delete) = false; %�޳�Q
        Q1(Delete) = true;  %����Q'
    end
    AC = Population(Choose);
    Rank = Rank(Choose);
    fS = fS(Choose);
end