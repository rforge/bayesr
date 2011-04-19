clear
infile intnr x pmean pqu2p5 pqu10 pmed pqu90 pqu97p5 pcat95 pcat80 using /tmp/Rtmpp1FwKj/bayesx/model_f_x_pspline.res
drop in 1
graph twoway rarea pqu2p5 pqu97p5 x, bcolor(gs13) || rarea pqu10 pqu90 x , bcolor(gs10) || /*
 */ scatter pmean x, c(l) m(i) clpattern(l) clcolor(gs0) /* 
 */ ytitle("Effect of x") xtitle("x") xlab(,grid) ylab(,grid) legend(off)
graph export /tmp/Rtmpp1FwKj/bayesx/model_f_x_pspline.eps, replace

