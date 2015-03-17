% usefile mcmc.prg

logopen using mcmc.prg.log

% gaussian test.
bayesreg b

dataset d 
d.infile using data.raw

b.outfile = mcmc

b.regress y = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3), family=gaussian iterations=12000 burnin=2000 step=10 setseed=1298345438 predict using d 

b.getsample

% poisson test.
%bayesreg b2
%b2.outfile = mcmc2

%b2.regress ycount = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3) + id(random), family=poisson iterations=12000 burnin=2000 step=10 setseed=1298345438 predict using d

%b2.getsample

logclose
