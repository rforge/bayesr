% usefile res_graphics.prg

dataset _dat
_dat.infile using res_f_x2_rw.res
graph _g
_g.plot x2 pmean pqu2p5 pqu10 pqu90 pqu97p5, title = "Effect of x2" xlab = x2 ylab = " " outfile = res_f_x2_rw.ps replace using _dat
drop _dat
drop _g

dataset _dat
_dat.infile using res_f_x1_pspline.res
graph _g
_g.plot x1 pmean pqu2p5 pqu10 pqu90 pqu97p5, title = "Effect of x1" xlab = x1 ylab = " " outfile = res_f_x1_pspline.ps replace using _dat
drop _dat
drop _g

dataset _dat
_dat.infile using res_f_district_spatial.res
map _map
_map.infile using input_filename
graph _g
_g.drawmap pmean district, map = _map color outfile = res_f_district_spatial_pmean.ps replace using _dat
_g.drawmap pcat95 district, map = _map nolegend pcat outfile = res_f_district_spatial_pcat95.ps replace using _dat
_g.drawmap pcat80 district, map = _map nolegend pcat outfile = res_f_district_spatial_pcat80.ps replace using _dat
drop _dat
drop _g
drop _map
