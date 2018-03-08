directoryname  = '2015-11-09/Binomial_lna_traj/';
experimentsprefix      = 'ParSet_traj_';
experimentsnumbers = 82:102;

K = 0.6:0.05:0.95;
binsnumbers = 25:25:300;

expnames = cell(max(size(experimentsnumbers)),1);
for i = 1:max(size(experimentsnumbers))
    expnames{i} = [experimentsprefix, num2str(experimentsnumbers(i))];
end

for expi = 1:max(size(expnames))
    expname = expnames{expi};
    directoryinput   = ['Input/',    directoryname];
    directoryoutput = ['Output/', directoryname];
    dirinput = [directoryinput,  expname, '.txt'];
    diroutput = [directoryoutput, expname, '/'];    
    directoryoutputresults = [diroutput, 'results/'];
    

     if ~exist([directoryoutputresults, 'B.csv'])
        expname = expnames{expi};
        ArimotoBlahutRun;
        ArimotoBlahutAnaliseResults;
     end
end