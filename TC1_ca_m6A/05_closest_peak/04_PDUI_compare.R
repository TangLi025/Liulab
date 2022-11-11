p0_pAPA_200 <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]

p5_pAPA_200 <- p5_pAPA_200[p5_pAPA_200$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]
p10_pAPA_200 <- p10_pAPA_200[p10_pAPA_200$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]
rp2_pAPA_200 <- rp2_pAPA_200[rp2_pAPA_200$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]

ggplot()+
  geom_density(data=p0_pAPA_200,mapping=aes(x=distance),color="red")+
  geom_density(data=p10_pAPA_200,mapping=aes(x=distance),color="green")

p0_pAPA_200_p10short <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="UP",]$Gene,1,18),] 
p0_pAPA_200_p10long <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="DOWN",]$Gene,1,18),]
p0_pAPA_200_p10nc <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="NC",]$Gene,1,18),] 

p10_pAPA_200_p10short <- p10_pAPA_200[p10_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="UP",]$Gene,1,18),] 
p10_pAPA_200_p10long <- p10_pAPA_200[p10_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="DOWN",]$Gene,1,18),] 
p10_pAPA_200_p10nc <- p10_pAPA_200[p10_pAPA_200$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="NC",]$Gene,1,18),] 


color=c("#fbefc4","#52696f","#e57b7f","#9e3150","#88baa4")
ggplot()+
  geom_density(data=p0_pAPA_200_p10long,mapping=aes(x=distance,color="p10long_p0_peak"))+
  geom_density(data=p0_pAPA_200_p10short,mapping=aes(x=distance,color="p10short_p0_peak"))+
  geom_density(data=p10_pAPA_200_p10long,mapping=aes(x=distance,color="p10long_p10_peak"))+
  geom_density(data=p10_pAPA_200_p10short,mapping=aes(x=distance,color="p10short_p10_peak"))+
  scale_color_manual(name="class",values = c("p10long_p0_peak"=color[5],
                                             "p10short_p0_peak"=color[2],
                                             "p10long_p10_peak"=color[3],
                                             "p10short_p10_peak"=color[4]))

ggplot()+
  geom_density(data=p10_pAPA_200_p10nc,mapping=aes(x=distance,color="p10nc_p10_peak"))+
  geom_density(data=p10_pAPA_200_p10long,mapping=aes(x=distance,color="p10long_p10_peak"))+
  geom_density(data=p10_pAPA_200_p10short,mapping=aes(x=distance,color="p10short_p10_peak"))+
  scale_color_manual(name="class",values = c("p10nc_p10_peak"=color[5],
                                             "p10long_p10_peak"=color[3],
                                             "p10short_p10_peak"=color[4]))

ggplot()+
  geom_density(data=p0_pAPA_200_p10nc,mapping=aes(x=distance,color="p10nc_p0_peak"))+
  geom_density(data=p0_pAPA_200_p10long,mapping=aes(x=distance,color="p10long_p0_peak"))+
  geom_density(data=p0_pAPA_200_p10short,mapping=aes(x=distance,color="p10short_p0_peak"))+
  scale_color_manual(name="class",values = c("p10nc_p0_peak"=color[5],
                                             "p10long_p0_peak"=color[3],
                                             "p10short_p0_peak"=color[4]))


ggplot()+
  geom_density(data=p10_pAPA_200_p10short,mapping=aes(x=distance),color="green")+
  geom_density(data=p0_pAPA_200_p10short,mapping=aes(x=distance),color="red")


#################### p0p5 ########################
ggplot()+
  geom_density(data=p0_pAPA_200,mapping=aes(x=distance),color="red")+
  geom_density(data=p5_pAPA_200,mapping=aes(x=distance),color="green")

p0_pAPA_200_p5short <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="UP",]$Gene,1,18),] 
p0_pAPA_200_p5long <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="DOWN",]$Gene,1,18),]
p0_pAPA_200_p5nc <- p0_pAPA_200[p0_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="NC",]$Gene,1,18),] 

p5_pAPA_200_p5short <- p5_pAPA_200[p5_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="UP",]$Gene,1,18),] 
p5_pAPA_200_p5long <- p5_pAPA_200[p5_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="DOWN",]$Gene,1,18),] 
p5_pAPA_200_p5nc <- p5_pAPA_200[p5_pAPA_200$name1 %in% substr(dapars_caRNA_p0p5[dapars_caRNA_p0p5$filter=="NC",]$Gene,1,18),] 


color=c("#fbefc4","#52696f","#e57b7f","#9e3150","#88baa4")
ggplot()+
  geom_density(data=p0_pAPA_200_p5long,mapping=aes(x=distance,color="p5long_p0_peak"))+
  geom_density(data=p0_pAPA_200_p5short,mapping=aes(x=distance,color="p5short_p0_peak"))+
  geom_density(data=p5_pAPA_200_p5long,mapping=aes(x=distance,color="p5long_p5_peak"))+
  geom_density(data=p5_pAPA_200_p5short,mapping=aes(x=distance,color="p5short_p5_peak"))+
  scale_color_manual(name="class",values = c("p5long_p0_peak"=color[5],
                                             "p5short_p0_peak"=color[2],
                                             "p5long_p5_peak"=color[3],
                                             "p5short_p5_peak"=color[4]))

ggplot()+
  geom_density(data=p5_pAPA_200_p5nc,mapping=aes(x=distance,color="p5nc_p5_peak"))+
  geom_density(data=p5_pAPA_200_p5long,mapping=aes(x=distance,color="p5long_p5_peak"))+
  geom_density(data=p5_pAPA_200_p5short,mapping=aes(x=distance,color="p5short_p5_peak"))+
  scale_color_manual(name="class",values = c("p5nc_p5_peak"=color[5],
                                             "p5long_p5_peak"=color[3],
                                             "p5short_p5_peak"=color[4]))

ggplot()+
  geom_density(data=p0_pAPA_200_p5nc,mapping=aes(x=distance,color="p5nc_p0_peak"))+
  geom_density(data=p0_pAPA_200_p5long,mapping=aes(x=distance,color="p5long_p0_peak"))+
  geom_density(data=p0_pAPA_200_p5short,mapping=aes(x=distance,color="p5short_p0_peak"))+
  scale_color_manual(name="class",values = c("p5nc_p0_peak"=color[5],
                                             "p5long_p0_peak"=color[3],
                                             "p5short_p0_peak"=color[4]))


ggplot()+
  geom_density(data=p5_pAPA_200_p5short,mapping=aes(x=distance),color="green")+
  geom_density(data=p0_pAPA_200_p5short,mapping=aes(x=distance),color="red")


######################### upstream500 #########
p0_pAPA_500 <- p0_pAPA_500[p0_pAPA_500$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]

p5_pAPA_500 <- p5_pAPA_500[p5_pAPA_500$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]
p10_pAPA_500 <- p10_pAPA_500[p10_pAPA_500$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]
rp2_pAPA_500 <- rp2_pAPA_500[rp2_pAPA_500$name1 %in% substr(dapars2_mean_omitNA$Gene,1,18),]

ggplot()+
  geom_density(data=p0_pAPA_500,mapping=aes(x=distance),color="red")+
  geom_density(data=p10_pAPA_500,mapping=aes(x=distance),color="green")

p0_pAPA_500_p10short <- p0_pAPA_500[p0_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="UP",]$Gene,1,18),] 
p0_pAPA_500_p10long <- p0_pAPA_500[p0_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="DOWN",]$Gene,1,18),]
p0_pAPA_500_p10nc <- p0_pAPA_500[p0_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="NC",]$Gene,1,18),] 

p10_pAPA_500_p10short <- p10_pAPA_500[p10_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="UP",]$Gene,1,18),] 
p10_pAPA_500_p10long <- p10_pAPA_500[p10_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="DOWN",]$Gene,1,18),] 
p10_pAPA_500_p10nc <- p10_pAPA_500[p10_pAPA_500$name1 %in% substr(dapars_caRNA_p0p10[dapars_caRNA_p0p10$filter=="NC",]$Gene,1,18),] 


color=c("#fbefc4","#52696f","#e57b7f","#9e3150","#88baa4")
ggplot()+
  geom_density(data=p0_pAPA_500_p10long,mapping=aes(x=distance,color="p10long_p0_peak"))+
  geom_density(data=p0_pAPA_500_p10short,mapping=aes(x=distance,color="p10short_p0_peak"))+
  geom_density(data=p10_pAPA_500_p10long,mapping=aes(x=distance,color="p10long_p10_peak"))+
  geom_density(data=p10_pAPA_500_p10short,mapping=aes(x=distance,color="p10short_p10_peak"))+
  scale_color_manual(name="class",values = c("p10long_p0_peak"=color[5],
                                             "p10short_p0_peak"=color[2],
                                             "p10long_p10_peak"=color[3],
                                             "p10short_p10_peak"=color[4]))

ggplot()+
  geom_density(data=p10_pAPA_500_p10nc,mapping=aes(x=distance,color="p10nc_p10_peak"))+
  geom_density(data=p10_pAPA_500_p10long,mapping=aes(x=distance,color="p10long_p10_peak"))+
  geom_density(data=p10_pAPA_500_p10short,mapping=aes(x=distance,color="p10short_p10_peak"))+
  scale_color_manual(name="class",values = c("p10nc_p10_peak"=color[5],
                                             "p10long_p10_peak"=color[3],
                                             "p10short_p10_peak"=color[4]))

ggplot()+
  geom_density(data=p0_pAPA_500_p10nc,mapping=aes(x=distance,color="p10nc_p0_peak"))+
  geom_density(data=p0_pAPA_500_p10long,mapping=aes(x=distance,color="p10long_p0_peak"))+
  geom_density(data=p0_pAPA_500_p10short,mapping=aes(x=distance,color="p10short_p0_peak"))+
  scale_color_manual(name="class",values = c("p10nc_p0_peak"=color[5],
                                             "p10long_p0_peak"=color[3],
                                             "p10short_p0_peak"=color[4]))


ggplot()+
  geom_density(data=p10_pAPA_500_p10short,mapping=aes(x=distance),color="green")+
  geom_density(data=p0_pAPA_500_p10short,mapping=aes(x=distance),color="red")
