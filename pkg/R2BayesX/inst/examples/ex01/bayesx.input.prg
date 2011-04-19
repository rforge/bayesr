% usefile /tmp/Rtmpp1FwKj/bayesx/bayesx.input.prg
dataset d
d.infile using /tmp/Rtmpp1FwKj/bayesx/bayesx.data.raw
bayesreg b
b.outfile = /tmp/Rtmpp1FwKj/bayesx/model
b.regress y = x(psplinerw2,nrknots=6,degree=3) + z*w(pspline2dimrw2,nrknots=9,degree=3) + fac2 + fac3 + fac4 + fac5 + fac6 + fac7 + fac8 + fac9 + fac10, family=gaussian iterations=1200 burnin=200 step=1 maxint=150 aresp=1 bresp=0.005 predict using d 
b.getsample
