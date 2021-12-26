function param = defaultParam()
    param = struct;
    param.method = 1;
    param.globalOptMethod = 0;
    param.maxIter = 2000;
    param.ftol = 1e-6;
    param.sampleNumPerIter = 50;
    param.gamma  = 0.2;
end

