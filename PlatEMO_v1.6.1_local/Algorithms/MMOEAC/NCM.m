function C=NCM(X,r)
    [n,dim]=size(X);
    S=(1:n);
    D=zeros(n,n);
    for i=1:dim
        D=D+(abs(X(:,i)-repmat(X(:,i)',n,1))<=r(i));
    end

    K=0;
    while size(S)>0 %当U不为空时
        K=K+1;
        Q=[];
        C(K).p=[];
        s=randperm(length(S));
        x=S(s(1)); %从U随机选择一个解
        Q=[Q x];%初始化Q
        C(K).p=[C(K).p x]; %初始化Ck
        while size(Q)>0 %当Q不为空时
            ss=randperm(length(Q));
            y=Q(ss(1));%随机选择Q中的解
            B=find(D(y,:)==dim);%取其直接相关解B(x)
            T=setdiff(B,C(K).p);%把B(x)中与Ck相同的解去除
            Q=[Q T];
            C(K).p=[C(K).p T];%将与直接相关解置入Q和聚类Ck中
            Q(ss(1))=[];%从Q中删除解x
        end
        S=setdiff(S,C(K).p); %把聚类Ck中的解从U中删除
    end
end