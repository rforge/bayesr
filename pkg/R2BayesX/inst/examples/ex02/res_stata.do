clear
infile intnr x2 pmean pqu2p5 pqu10 pmed pqu90 pqu97p5 pcat95 pcat80 using res_f_x2_rw.res
drop in 1
graph twoway rarea pqu2p5 pqu97p5 x2, bcolor(gs13) || rarea pqu10 pqu90 x2 , bcolor(gs10) || /*
 */ scatter pmean x2, c(l) m(i) clpattern(l) clcolor(gs0) /* 
 */ ytitle("Effect of x2") xtitle("x2") xlab(,grid) ylab(,grid) legend(off)
graph export res_f_x2_rw.eps, replace

clear
infile intnr x1 pmean pqu2p5 pqu10 pmed pqu90 pqu97p5 pcat95 pcat80 using res_f_x1_pspline.res
drop in 1
graph twoway rarea pqu2p5 pqu97p5 x1, bcolor(gs13) || rarea pqu10 pqu90 x1 , bcolor(gs10) || /*
 */ scatter pmean x1, c(l) m(i) clpattern(l) clcolor(gs0) /* 
 */ ytitle("Effect of x1") xtitle("x1") xlab(,grid) ylab(,grid) legend(off)
graph export res_f_x1_pspline.eps, replace

