source("low_rhc_cont_sim.r")
source("low_rhc_bin_sim.r")
source("low_studentperf_cont_sim.r")
source("low_studentperf2_cont_sim.r")
source("low_cervicalcancer_bin_sim.r")

set.seed(10000)
ix<-sample(1:2000,2000,F)

for(i in seq(0,1980,by=20))
 {
low_rhc_cont_sim(ix[1+i],ix[2+i],ix[3+i],ix[4+i])
low_studentperf_cont_sim(ix[5+i],ix[6+i],ix[7+i],ix[8+i])
low_studentperf2_cont_sim(ix[9+i],ix[10+i],ix[11+i],ix[12+i])
low_rhc_bin_sim(ix[13+i],ix[14+i],ix[15+i],ix[16+i])
low_cervicalcancer_bin_sim(ix[17+i],ix[18+i],ix[19+i],ix[20+i])
print(i/20)
 }