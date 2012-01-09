% usefile reml.prg

logopen using reml.prg.log

remlreg b

dataset d 
d.infile using data.raw

b.outfile = reml

b.regress y = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3), family=gaussian eps=1e-05 lowerlim=0.001 maxit=400 maxchange=1e+06 using d 

logclose 
