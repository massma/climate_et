## Meet with Pierre ##

After meeting with Pierre he said a couple things would be useful. First, looking at when soil moisture variability is important vs when VPD variability is improtant. The approach to this would be somehow normalizing by soil moisture variability and VPD varibility. (right, because even if dET/dVPD is high, if VPD never varies at a site it doesn't matter.). What is difficult though is that soil moisture might only matter at certain time scales or lags. Julia is working on this now so her results might help me understand this. So maybe compare both dET/dVPD normalized by std(VPD), and likewise for soil moisture. One idea to get an analytical formulation for soil moisture is to fit an empiracle funtion (many people have used this in the past, pierre thought the weibull fit was one that was used), and then calculate analytical derivatives off of the empiracle fits.

A second paper building off this project would be to then create maps where soil moisture variability dominates vs atmospheric demand. For example, there is quite a bit fo certainty that soil in the horn of Africa and the Mediteranean will become much drier.  So there soil moisture variability will dominate, but atmospheric demand will probalby dominate anywhere there are crops (and probably everywehre else).

Also, looking at dET/dLW would be interesting as much of the change in radiation with GW should be in the LW.

However, Pierre mentioned he thinks there is a LW-SW feedback. This makes sense: for example LW increases surface T, which decreases evaporative fraction, which decreases clouds, and increases SW and PAR. Looking at this could be an interesting futre project.

Also, Pierre said to just use z0 = 0.1*h for surface roughness. The LAI term is negligable but would make the analytical formulation much more difficult.