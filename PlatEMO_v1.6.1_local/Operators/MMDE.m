function Offspring = MMDE(Global,Parent)
% <operator> <real>
% Differental evolution and polynomial mutation
% CR   ---   1 --- Parameter CR in differental evolution
% F    --- 0.5 --- Parameter F in differental evolution
% proM ---   1 --- The expectation of number of bits doing mutation 
% disM ---  20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [CR,F,proM,disM] = Global.ParameterSet(1,0.5,1,20);
    Parent    = Parent([1:end,1:ceil(end/5)*5-end]);
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);

    %% Differental evolution  
    % 输入规模为5N
    Parent1Dec   = ParentDec(1:N/5,:);
    Parent2Dec   = ParentDec(N/5+1:N/5*2,:);
    Parent3Dec   = ParentDec(N/5*2+1:N/5*3,:);
    Parent4Dec   = ParentDec(N/5*3+1:N/5*4,:);
    Parent5Dec   = ParentDec(N/5*4+1:end,:);
    OffspringDec1 = Parent1Dec;
    Site = rand(N/5,D) < CR;
    OffspringDec1(Site) = Parent1Dec(Site) + F*(Parent2Dec(Site)-Parent3Dec(Site)) + F*(Parent4Dec(Site)-Parent5Dec(Site));
%     OffspringDec2 = Parent2Dec;
%     Site = rand(N/5,D) < CR;
%     OffspringDec2(Site) = Parent2Dec(Site) + F*(Parent1Dec(Site)-Parent3Dec(Site)) + F*(Parent4Dec(Site)-Parent5Dec(Site));
%     OffspringDec3 = Parent3Dec;
%     Site = rand(N/5,D) < CR;
%     OffspringDec3(Site) = Parent3Dec(Site) + F*(Parent1Dec(Site)-Parent2Dec(Site)) + F*(Parent4Dec(Site)-Parent5Dec(Site));
%     OffspringDec4 = Parent4Dec;
%     Site = rand(N/5,D) < CR;
%     OffspringDec4(Site) = Parent4Dec(Site) + F*(Parent1Dec(Site)-Parent2Dec(Site)) + F*(Parent3Dec(Site)-Parent5Dec(Site));
%     OffspringDec5 = Parent5Dec;
%     Site = rand(N/5,D) < CR;
%     OffspringDec5(Site) = Parent5Dec(Site) + F*(Parent1Dec(Site)-Parent2Dec(Site)) + F*(Parent3Dec(Site)-Parent4Dec(Site));
%     
   
    %% Polynomial mutation
     Lower = repmat(Global.lower,N,1);
     Upper = repmat(Global.upper,N,1);
     Site  = rand(N/5,D) < proM/D;
     mu    = rand(N/5,D);
     temp  = Site & mu<=0.5;
     OffspringDec1(temp) = OffspringDec1(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
                         (1-(OffspringDec1(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
     temp = Site & mu>0.5; 
     OffspringDec1(temp) = OffspringDec1(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
                         (1-(Upper(temp)-OffspringDec1(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
     
%      Site  = rand(N/5,D) < proM/D;
%      mu    = rand(N/5,D);
%      temp  = Site & mu<=0.5;
%      OffspringDec2(temp) = OffspringDec2(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                          (1-(OffspringDec2(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%      temp = Site & mu>0.5; 
%      OffspringDec2(temp) = OffspringDec2(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                          (1-(Upper(temp)-OffspringDec2(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%      
%      Site  = rand(N/5,D) < proM/D;
%      mu    = rand(N/5,D);
%      temp  = Site & mu<=0.5;
%      OffspringDec3(temp) = OffspringDec3(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                          (1-(OffspringDec3(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%      temp = Site & mu>0.5; 
%      OffspringDec3(temp) = OffspringDec3(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                          (1-(Upper(temp)-OffspringDec3(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%                      
%      Site  = rand(N/5,D) < proM/D;
%      mu    = rand(N/5,D);
%      temp  = Site & mu<=0.5;
%      OffspringDec4(temp) = OffspringDec4(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                          (1-(OffspringDec4(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%      temp = Site & mu>0.5; 
%      OffspringDec4(temp) = OffspringDec4(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                          (1-(Upper(temp)-OffspringDec4(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));                
%                      
%      Site  = rand(N/5,D) < proM/D;
%      mu    = rand(N/5,D);
%      temp  = Site & mu<=0.5;
%      OffspringDec5(temp) = OffspringDec5(temp)+(Upper(temp)-Lower(temp)).*((2.*mu(temp)+(1-2.*mu(temp)).*...
%                          (1-(OffspringDec5(temp)-Lower(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1))-1);
%      temp = Site & mu>0.5; 
%      OffspringDec5(temp) = OffspringDec5(temp)+(Upper(temp)-Lower(temp)).*(1-(2.*(1-mu(temp))+2.*(mu(temp)-0.5).*...
%                          (1-(Upper(temp)-OffspringDec5(temp))./(Upper(temp)-Lower(temp))).^(disM+1)).^(1/(disM+1)));
%      OffspringDec = [OffspringDec1;OffspringDec2;OffspringDec3;OffspringDec4;OffspringDec5];
%      OffspringDec = [OffspringDec1;OffspringDec2];
      OffspringDec = real(OffspringDec1);
     for i = 1:size(OffspringDec,1)
         if sum(sum(OffspringDec(i,:) >= Global.upper)) > 0 || sum(sum(OffspringDec(i,:) <= Global.lower)) > 0 
             %Site = rand(N/5,D) < 1;
             d = i;
             while d > N/5
                 d = d - N/5;
             end
             OffspringDec(i) = Parent1Dec(d) - F*(Parent2Dec(d)-Parent3Dec(d)) - F*(Parent4Dec(d)-Parent5Dec(d));
           if sum(sum(OffspringDec(i,:) >= Global.upper)) > 0 || sum(sum(OffspringDec(i,:) <= Global.lower)) > 0 
                for j = 1:length(OffspringDec(i,:))
                    if OffspringDec(i,j) >= Global.upper(j)
                        OffspringDec(i,j) = Global.lower(j) + (OffspringDec(i,j) - Global.upper(j));
                    end
                    if OffspringDec(i,j) <= Global.lower(j)
                        OffspringDec(i,j) = Global.upper(j) - (Global.lower(j) - OffspringDec(i,j));
                    end
                end
            end
         end
      end
%     for i = 1:length(OffspringDec)
%         temp1 = INDIVIDUAL(OffspringDec(i,:));
%         temp2 = INDIVIDUAL(Parent1Dec(i,:));
%         [FrontNo,MaxFNo] = NDSort([temp1.objs;temp2.objs],inf);
%         choose = FrontNo == MaxFNo; %(0,1)表示父代支配新解
%         if choose(1) > choose(2)
%             OffspringDec(i,:) = Parent1Dec(i,:);
%         end
%     end
    
%     
    Offspring = INDIVIDUAL(OffspringDec);
end