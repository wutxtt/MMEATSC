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
   
    Xic  = false(1,Global.D); %���������Ա���
    Xre  = false(1,Global.D); %
    
    %% Find convergence-related variable
    for i = 1 : Global.D
        x      = rand(1,Global.D).*(Global.upper-Global.lower) + Global.lower;
        S      = repmat(x,NCA,1); %NCA*D����
        S(:,i) = ((1:NCA)'-1+rand(NCA,1))/NCA*(Global.upper(i)-Global.lower(i)) + Global.lower(i);
        [~,S]  = Global.problem('value',Global,S); %���ʼ����ı��i������������Ŀ��ֵ
        S = unique(S,'rows'); % delete the duplicate ɾ���ظ��������˳�򷵻�S���ظ���Ŀ��ֵ
        [~,MaxFNo] = NDSort(S,inf);
        if MaxFNo == size(S,1) %��֧������ı���
            Xic(i) = true; %��Ϊ�������������
        else
            Xre(i) = true; %����Ϊ����Ա���
        end
    end
    
    %% Interdependence �����
    % Generate the initial population        
    PopDec = rand(Global.N,Global.D);
    PopDec = PopDec.*repmat(Global.upper-Global.lower,Global.N,1) + repmat(Global.lower,Global.N,1); 
    %�������һ�����壬��������Ŀ��ֵ
    [~,PopObj]= Global.problem('value',Global,PopDec);
    % Interdependence analysis
    interaction = false(Global.D); %��ʼ�������߼�����
    interaction(logical(eye(Global.D))) = true; %��λ����ת��Ϊ�߼�ֵ��ÿ�����������������
    for i = 1 : Global.D-1
        for j = i+1 : Global.D
            for time2try = 1 : NIA
                % Detect whether the i-th and j-th decision variables are interacting
                x    = randi(Global.N);%1~N֮����������
                a2   = rand*(Global.upper(i)-Global.lower(i)) + Global.lower(i); %���ȡ��i������
                b2   = rand*(Global.upper(j)-Global.lower(j)) + Global.lower(j); %���ȡ��j������
                Decs = repmat(PopDec(x,:),3,1); %���ȡһ������ľ��߿ռ���������������
                Decs(1,i) = a2;
                Decs(2,j) = b2;
                Decs(3,[i,j]) = [a2,b2]; %�����i������j���������������ֱ��滻���Ƚ�
                [~,F]= Global.problem('value',Global,Decs); %�ֱ�õ�����Ŀ��ֵ
                delta1 = F(1,:) - PopObj(x,:); %jֵ������ı��i�����߱������ÿ��Ŀ��ֵ�ı仯��
                delta2 = F(3,:) - F(2,:); %����һ����ͬjֵ�ı��i�����߱���
                interaction(i,j) = interaction(i,j) | any(delta1.*delta2<0); 
                %�仯�������Ŀ��ֵ��������Ŀ��ֵ�仯�̶��෴��Ϊ���
                interaction(j,i) = interaction(i,j);                
            end
        end
    end
    
    %% Group based on Interdependence
    while sum(sum(interaction(Xic,Xre)))
        for i = find(Xic==1) %�ҵ�һ�����߱����Ƕ�����
            fprintf('i=%d\n',i);
            if sum(interaction(i,Xre)) %������ĳ���������
                Xic(i) = false; %������������
                Xre(i) = true; %��������������洢����
            end
        end      
    end

end