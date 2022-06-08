function Offspring = Gaussian(Global,Parent)
% <operator> <real>
% Simulated binary crossover and polynomial mutation
% proC ---  1 --- The probability of doing crossover
% disC --- 20 --- The distribution index of simulated binary crossover
% proM ---  1 --- The expectation of number of bits doing mutation 
% disM --- 20 --- The distribution index of polynomial mutation

%--------------------------------------------------------------------------
% Copyright (c) 2016-2017 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB Platform
% for Evolutionary Multi-Objective Optimization [Educational Forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    [u, sigma] = Global.ParameterSet(0,2); %高斯分布参数 
    
    ParentDec = Parent.decs;
    [N,D]     = size(ParentDec);
    gaussian  = sigma * randn(1,D) + u;
    OffspringDec = ParentDec + gaussian;

    Offspring = INDIVIDUAL(OffspringDec);
end