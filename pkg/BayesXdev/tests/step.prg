% usefile step.prg

logopen using step.prg.log

stepwisereg b

dataset d 
d.infile using data.raw

b.outfile = step

b.regress y = x1(psplinerw2,nrknots=20,degree=3) + x2(psplinerw2,nrknots=20,degree=3) + x3(psplinerw2,nrknots=20,degree=3) + x4(psplinerw2,nrknots=20,degree=3), family=gaussian iterations=10000 burnin=2000 step=10 setseed=1298345438 algorithm=cdescent1 CI=MCMCselect predict using d 

logclose 
