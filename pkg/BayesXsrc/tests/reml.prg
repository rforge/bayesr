% usefile reml.prg

logopen using reml.prg.log

% gaussian test
remlreg b

dataset d 
d.infile using data.raw

b.outfile = reml

b.regress y = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3), family=gaussian eps=1e-05 lowerlim=0.001 maxit=400 maxchange=1e+06 using d 

% poisson test.
%remlreg b2
%b2.outfile = reml2

%b2.regress ycount = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3) + id(random), family=poisson eps=1e-05 lowerlim=0.001 maxit=400 maxchange=1e+06 using d 

%b2.getsample

logclose
