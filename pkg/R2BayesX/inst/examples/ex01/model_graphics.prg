% usefile /tmp/Rtmpp1FwKj/bayesx/model_graphics.prg

dataset _dat
_dat.infile using /tmp/Rtmpp1FwKj/bayesx/model_f_x_pspline.res
graph _g
_g.plot x pmean pqu2p5 pqu10 pqu90 pqu97p5, title = "Effect of x" xlab = x ylab = " " outfile = /tmp/Rtmpp1FwKj/bayesx/model_f_x_pspline.ps replace using _dat
drop _dat
drop _g
