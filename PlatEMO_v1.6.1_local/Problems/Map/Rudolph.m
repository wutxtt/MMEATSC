function varargout = Rudolph(Operation,Global,input)
% <problem> <MMPPO>
% Multi-modal Multi-objective test Function
% operator --- EAreal

%--------------------------------------------------------------------------
% Copyright 2017-2018 Yiping Liu
% This is the code of MMF used in "Yiping Liu, Gary G. Yen, 
% and Dunwei Gong, A Multi-Modal Multi-Objective Evolutionary Algorithm 
% Using Two-Archive and Recombination Strategies, IEEE Transactions on 
% Evolutionary Computation, 2018, Early Access".
% Please contact {yiping0liu@gmail.com} if you have any problem.
%--------------------------------------------------------------------------
% MMF is proposed in " Caitong Yue, Boyang Qu, and Jing Liang, 
% A Multi-objective Particle Swarm Optimizer Using Ring Topology for 
% Solving Multimodal Multi-objective Problems, IEEE Transactions on 
% Evolutionary Computation, 2017, Early Access".
%--------------------------------------------------------------------------
% This code uses PlatEMO published in "Ye Tian, Ran Cheng, Xingyi Zhang, 
% and Yaochu Jin, PlatEMO: A MATLAB Platform for Evolutionary 
% Multi-Objective Optimization [Educational Forum], IEEE Computational 
% Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    switch Operation
        case 'init'
            Global.M        = 2;
            Global.D        = 2;
            Global.lower    = [-8,-8];
            Global.upper    = [8,8];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            a = 0.5;
            b = 5;
            c = 5;
            t1 = sign(PopDec(:,1)) .* min(ceil((abs(PopDec(:,1))-a-c/2)/(2*a+c)),1);
            t2 = sign(PopDec(:,2)) .* min(ceil((abs(PopDec(:,2))-b/2)/b),1);
            f = t1==t2;
            ft = ~f==(t1*100);
            ft = ~ft*0.1;
            PopObj = NaN(N,Global.M);
            PopObj(:,1) = (PopDec(:,1) - t1.*((c+2*a)+a)).^2 + (PopDec(:,2) - t2.*b).^2 + ft;
            PopObj(:,2) = (PopDec(:,1) - t1.*((c+2*a)-a)).^2 + (PopDec(:,2) - t2.*b).^2 + ft;         
            
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = [];
            varargout = {f};        
    end
end