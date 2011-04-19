% usefile c:\arbeit\packages\bayesx\examples\examples.prg

dataset d
bayesreg b
map m

% nonparametric

d.infile using c:\arbeit\packages\bayesx\inst\examples\nonparametric.raw

b.outfile = c:\arbeit\packages\bayesx\inst\examples\nonparametric
b.regress y = x(psplinerw2), family=gaussian iterations=12000 burnin=2000 step=10 using d
b.getsample

b.outfile = c:\arbeit\packages\bayesx\inst\examples\nonparametric2
b.regress ytime = time(psplinerw2), family=gaussian iterations=12000 burnin=2000 step=10 using d

% spatial

d.infile using c:\arbeit\packages\bayesx\inst\examples\spatial.raw
m.infile using c:\arbeit\packages\bayesx\inst\examples\germany.bnd

b.outfile = c:\arbeit\packages\bayesx\inst\examples\spatial
b.regress y = regions(spatial, map=m), family=gaussian iterations=12000 burnin=2000 step=10 using d


% surface

d.infile using c:\arbeit\packages\bayesx\inst\examples\surface.raw

b.outfile = c:\arbeit\packages\bayesx\inst\examples\surface
b.regress y = x1*x2(pspline2dimrw1, nrknots=12), family=gaussian iterations=12000 burnin=2000 step=10 using d

