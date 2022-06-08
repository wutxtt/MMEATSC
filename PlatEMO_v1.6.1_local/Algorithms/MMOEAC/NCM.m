function C=NCM(X,r)
    [n,dim]=size(X);
    S=(1:n);
    D=zeros(n,n);
    for i=1:dim
        D=D+(abs(X(:,i)-repmat(X(:,i)',n,1))<=r(i));
    end

    K=0;
    while size(S)>0 %��U��Ϊ��ʱ
        K=K+1;
        Q=[];
        C(K).p=[];
        s=randperm(length(S));
        x=S(s(1)); %��U���ѡ��һ����
        Q=[Q x];%��ʼ��Q
        C(K).p=[C(K).p x]; %��ʼ��Ck
        while size(Q)>0 %��Q��Ϊ��ʱ
            ss=randperm(length(Q));
            y=Q(ss(1));%���ѡ��Q�еĽ�
            B=find(D(y,:)==dim);%ȡ��ֱ����ؽ�B(x)
            T=setdiff(B,C(K).p);%��B(x)����Ck��ͬ�Ľ�ȥ��
            Q=[Q T];
            C(K).p=[C(K).p T];%����ֱ����ؽ�����Q�;���Ck��
            Q(ss(1))=[];%��Q��ɾ����x
        end
        S=setdiff(S,C(K).p); %�Ѿ���Ck�еĽ��U��ɾ��
    end
end