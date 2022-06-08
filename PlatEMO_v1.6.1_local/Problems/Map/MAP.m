function varargout = MAP(Operation,Global,input)
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
            Global.M        = 4;
            Global.D        = 2;
            Global.lower    = [0,0];
            Global.upper    = [100,100];
            Global.operator = @EAreal;

            PopDec    = rand(input,Global.D).*repmat(Global.upper-Global.lower,input,1)+repmat(Global.lower,input,1);
            varargout = {PopDec};
        case 'value'
            Point.Es=[3,37;42,96;45,60;50,25;83,72;98,38];
            Point.Js=[40,20;51,60;95,51];
            Point.Cs=[10,55;15,15;15,78;15,88;20,23;20,70;32,42;35,60;40,76;52,78;52,96;55,33;75,27];
            Point.Rs=[17.5,82.5;55.5,82.5;94.5,6.5];
%             Es=[3,37;42,96;45,60;50,25;83,72;98,38];
%             Js=[40,20;51,60;95,51];
%             Cs=[10,55;15,15;15,78;15,88;20,23;20,70;32,42;35,60;40,76;52,78;52,96;55,33;75,27];
%             Rs=[17.5,82.5;55.5,82.5;94.5,6.5];
%             Point = [Es; Js; Cs; Rs];
            PopDec = input;
            [N,~]  = size(PopDec);
            PopObj = zeros(N,Global.M);
            for i=1:N
                for p=1:size(Point.Es,1)
                    X = PopDec(i,:);
                    Y = Point.Es;
                    DEs = [];
                    for y=1:size(Y,1)
                        DEs = [DEs, pdist2(X, Y(y,:))];
%                         PopObj(i,1) = PopObj(i,1) + pdist2(X, Y(y,:));
                    end
                    PopObj(i,1) = min(DEs);
                end
                for p=1:size(Point.Js,1)
                    X = PopDec(i,:);
                    Y = Point.Js;
                    DJs = [];
                    for y=1:size(Y,1)
                        DJs = [DJs, pdist2(X, Y(y,:))];
%                         PopObj(i,1) = PopObj(i,1) + pdist2(X, Y(y,:));
                    end
                    PopObj(i,2) = min(DJs);
                end
                for p=1:size(Point.Cs,1)
                    X = PopDec(i,:);
                    Y = Point.Cs;
                    DCs = [];
                    for y=1:size(Y,1)
                        DCs = [DCs, pdist2(X, Y(y,:))];
%                         PopObj(i,1) = PopObj(i,1) + pdist2(X, Y(y,:));
                    end
                    PopObj(i,3) = min(DCs);
                end
                for p=1:size(Point.Rs,1)
                    X = PopDec(i,:);
                    Y = Point.Rs;
                    DRs = [];
                    for y=1:size(Y,1)
                        DRs = [DRs, pdist2(X, Y(y,:))];
%                         PopObj(i,1) = PopObj(i,1) + pdist2(X, Y(y,:));
                    end
                    PopObj(i,4) = min(DRs);
                end
            end
%             PopObj(:,1) = abs(PopDec(:,1)-2);
%             PopObj(:,2) = 1-sqrt(PopObj(:,1))+2*(PopDec(:,2)-sin(6*pi*PopObj(:,1)+pi)).^2;
%             
            PopCon = [];
    
            varargout = {input,PopObj,PopCon};
        case 'PF'
            f = [];
            varargout = {f};        
    end
end