% remove project directories to matlab search path

rmpath(genpath([pwd '/algorithm/']));
rmpath(genpath([pwd '/external/']));
rmpath(genpath([pwd '/test/']));
rmpath(genpath([pwd '/utility/']));
rmpath(pwd);
savepath