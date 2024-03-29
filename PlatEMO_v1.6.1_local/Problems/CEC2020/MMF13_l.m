function varargout = MMF13(Operation,Global,input)
% <problem> <MMF>
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
            Global.D        = 3;
            Global.lower    = [0.1 0.1 0.1];
            Global.upper    = [1.1 1.1 1.1];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = NaN(N,Global.M);
            PopObj(:,1) = PopDec(:,1);
            g=2-exp(-2*log10(2)*((PopDec(:,2)+sqrt(PopDec(:,3))-0.1)/0.8).^2).*(sin(2*pi*(PopDec(:,2)+sqrt(PopDec(:,3))))).^6;
            PopObj(:,2) = g./PopDec(:,1);
            
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f =[];
            varargout = {f};        
    end
end