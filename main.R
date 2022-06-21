#------------Construct models--------
source('funcitons.R')
library(spatstat)
library(igraph)
library(dplyr)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
#set the parameters
ground <- owin(xrange = c(0, 200), yrange = c(0,200)) #set the investigation area
max_crown <- 5
min_crown <- 2
lambda <- 0.015 #
radius <- 5  #Diffusion radius
HC_R<-4 #HC model
ave_offspring_per_cluster <- 3

#------------Monte-Carlo simulations--------
type='WCL';#'CS'/'CL'/'WCL'

n_simu=199
directory='./simu_data/'
set.seed(123)
ttype=data.frame()
#ttype=read.csv(paste0(directory,'/','all_index.csv'))
for (type in c('CS','CL','WCL')){
cat('=============!! running ',type,' !!==============\n')
dir.create(paste0(directory,type))
tmod=data.frame()
  for (mod in c('Tho','Mat','CSR','Str','HC')){
    cat('=============running ',mod,' ==============\n')
    indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                       mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
    for (j in 1:n_simu){
      make_mod(mod = mod)->mod1
      make_net(mod1,type = type)->tree_net
      write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
      net_par(tree_net,type = type)->net_indx
      indexs<-rbind(indexs,net_indx)
      if (j%%10==0)cat(j,'have done\n')
    }
    tmod=rbind(all_index,data.frame(mod=mod,num=1:n_simu,indexs))
  }
ttype=rbind(ttype,data.frame(net=type,tmod))
}
write.csv(ttype,file = paste0(directory,'/','all_index.csv'),row.names = F)


set.seed(123)
modspls=list()
for (mod in c('Tho','Mat','CSR','Str','HC')){
  make_mod(mod = mod)->mod1
  modspls[[paste0(mod,'1')]]=plot_mod(mod1)[[1]]+labs(subtitle = mod)+scale_size(range=c(0.5,2))+theme_pubr(legend = 'right')
  modspls[[paste0(mod,'2')]]=plot_mod(mod1)[[2]]+labs(subtitle = mod)+scale_size(range=c(0.5,2))+theme_pubr(legend = 'right')
  
}
cowplot::plot_grid(plotlist = modspls,align = 'hv',labels = 'auto',nrow = 5,label_size = 14)%>%
  ggsave('figures/mods.pdf',plot = .,width =9,height = 16 )

pdf('figures/fig2.pdf',width = 20,height = 12)
graphics::layout(matrix(1:15,3,5,byrow=F))
for (mod in c('Tho','Mat','CSR','Str','HC')){
  par(mar=c(0,0,1,0)) 
  make_mod(mod = mod)->mod1
  plot_net(mod1,type ='CS',main = mod,size_thr = 4.9)
  plot_net(mod1,type ='CL',main = mod,size_thr = 4.9)
  plot_net(mod1,type ='WCL',main = mod,size_thr = 4.9)
}
dev.off()

cowplot::plot_grid(plotlist = modspls,align = 'hv',labels = 'auto',nrow = 5,label_size = 14)%>%
  ggsave('figures/mods.pdf',plot = .,width =9,height = 16 )

#-------------Box_plot--------------------
library(dplyr)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(patchwork)
library(RColorBrewer)

all_index<-do.call(rbind,list(data.frame(net='CS',read.csv('simu_data/CS.csv')),
                              data.frame(net='CL',read.csv('simu_data/CL.csv')),
                              data.frame(net='WCL',read.csv('simu_data/WCL.csv'))))

#write.csv(all_index,'simu_data/all_index.csv',row.names = F)
#read.csv('simu_data/all_index.csv')->all_index

all_index$net<-factor(all_index$net,levels = c('CS','CL','WCL'))
all_index$mod<-factor(all_index$mod,levels = c('Tho','Mat','CSR','Str','HC'))
#all indexs in five models of three networks
melt(all_index,id.vars = c('net','mod','num'),variable.name = 'indexs')->long_index
if (F){
  lsd<-data.frame()
  for (net in levels(long_index$net)){ 
    for (indexs in levels(long_index$indexs)){
      long_index%>%filter(net=={!!net},indexs=={!!indexs})->tmpdf
      lsd<-rbind(lsd,data.frame(net=net,indexs=indexs,multitest(tmpdf$value,tmpdf$mod)))
    }
  }
  lsd<-mutate(lsd,net=factor(net,levels = levels(long_index$net)),indexs=factor(indexs,levels = levels(long_index$indexs)))
  lsd%>%group_by(net,indexs)%>%mutate(high=max(var))->lsd
  
  ap<-ggboxplot(long_index,x = 'mod',y='value',color='mod',add = c('jitter'),
            add.params = list(color='mod',width = 0.15,alpha=0.8,size=0.5),outlier.shape = NA)+
    scale_color_manual(values =brewer.pal(9,'Set1'))+ylab(label = NULL)+xlab(label = NULL)+
    geom_text(data = lsd,aes(x=variable,y=1.2*high,label=groups), inherit.aes = FALSE ,color='red',size=5)+
    facet_wrap(net~indexs,scales = 'free_y',nrow = 3)+
    theme_pubr(base_size = 14)+ theme(plot.margin=unit(rep(1,4),'lines'),legend.position="none")
  pdf('figures/all_box.pdf',35,15)
  ap
  dev.off()
}
select_nets=c('CS','CL','WCL')
select_inds=c('mdegree','clusteringC','edensity','aplength')
anodf<-list()
for (net in select_nets){
  for (indexs in select_inds){
    long_index%>%filter(net=={!!net},indexs=={!!indexs})->tmpdf
    aov(value~mod,data = tmpdf)->ano
    p=summary(ano)[[1]][1,5]
    names(p)<-'anova'
    TukeyHSD(ano, 'mod', p.adj = 'bonferroni')$mod->aa
    aa=c(p,aa[,4])
    data.frame(aa)->aa
    colnames(aa)<-paste(net,indexs,sep = '_')
    anodf[[paste(net,indexs,sep = '_')]]<-aa
  }
}
do.call(cbind,anodf)%>%signif(.,3)%>%write.csv(.,'qw.csv')

select_nets=c('CS','CL','WCL')
select_inds=c('mdegree','clusteringC','edensity','aplength')

boxpls<-list()
for (net in select_nets){
  for (indexs in select_inds){
    boxpls[[paste0(net,'_',indexs)]]<-box_plot(long_index,net,indexs)
  }
}

pdf('figures/fig3.pdf',16,12)
cowplot::plot_grid(plotlist = boxpls,nrow = 3,align ='hv',labels ='auto',label_size = 21)+
  theme(plot.margin=unit(c(3,0,0,5),'lines'))
dev.off()

#summary
all_index%>%group_by(net,mod)%>%summarise(min(nedges),max(nedges),)
long_index%>%group_by(indexs,net,mod)%>%summarise(Max = max(value),
                                                  Min=min(value),
                                                  Mean=mean(value),
                                                  SD=sd(value),
                                                  CV=sd(value)/mean(value))->stat_res
write.csv(stat_res,'simu_data/stat_res.csv',row.names = F)

all_index%>%group_by(net,mod)%>%summarise(nnodes=paste(round(min(nnodes),0),round(mean(nnodes),0),round(max(nnodes),0),sep = '-'),
                                          nedges=paste(round(min(nedges),0),round(mean(nedges),0),round(max(nedges),0),sep = '-'),
                                          k=paste(round(min(mdegree),2),round(mean(mdegree),2),round(max(mdegree),2),sep = '-'),
                                          C=paste(round(min(clusteringC),2),round(mean(clusteringC),2),round(max(clusteringC),2),sep = '-'),
                                          D=paste(round(min(edensity),4),round(mean(edensity),4),round(max(edensity),4),sep = '-'),
                                          L=paste(round(min(aplength),2),round(mean(aplength),2),round(max(aplength),2),sep = '-'))%>%
  write.csv(.,'simu_data/stat_res_2.csv',row.names = F)

#-------------Distribution--------------------
library(ggridges)
library(viridis)
library(hrbrthemes)


#distribution of degree or weighted edges
type='CL';#'CS'/'CL'/'WCL'
directory='./simu_data/distrbution/'
set.seed(123)
disofd=data.frame();disofe=data.frame()
for (mod in c('Tho','Mat','CSR','Str','HC')){
  cat('=============running ',mod,' ==============\n')
  make_mod(mod = mod)->mod1
  make_net(mod1,type = type)->tree_net
  
  if (!type%in%c('CS','CL','WCL'))stop('type is must be one of CS,CL,WCL')
  if(type%in%c('CS','CL')){
    p<-graph_from_data_frame(tree_net,directed = F)
    disofd=rbind(disofd,data.frame(mod=mod,degree=as.vector(degree(p))))
  }
  if(type=='WCL'){
    wdegree<-\(p){
      a=degree(p);b=E(p)$weight
      wd=c();index=1
      for (i in a/2){
        tmp=b[index:(index+i-1)]
        wd=c(wd,sum(tmp))
        index=index+i
      }
      return(wd) 
    }
    p<-graph_from_data_frame(tree_net)
    disofd=rbind(disofd,data.frame(mod=mod,degree=wdegree(p)))
    disofe=rbind(disofe,data.frame(mod=mod,weight=E(p)$weight))
  }
}
write.csv(disofd,file = paste0(directory,'/',type,'_disofd.csv'),row.names = F)
if(type=='WCL')write.csv(disofe,file = paste0(directory,'/',type,'_disofe.csv'),row.names = F)

if(T){
  dis1<-read.csv('./simu_data/distrbution/CS_disofd.csv')
  dis2<-read.csv('./simu_data/distrbution/CL_disofd.csv')
  dis3<-read.csv('./simu_data/distrbution/WCL_disofd.csv')
  dis4<-read.csv('./simu_data/distrbution/WCL_disofe.csv')
}


p1<-displot(dis1)+
  labs(title = expression('CS distribution of degrees'))
p2<-displot(dis2)+
  labs(title = expression('CL distribution of degrees'))

p3<-displot(dis3)+
  labs(title = expression('WCL distribution of degrees'))+xlab(label = 'wdegree')

p4<-displot(dis4)+
  labs(title = expression('WCL distribution of CI'[ij]))+xlab(label = expression('CI'[ij]))

pdf('figures/fig4.pdf',width = 8,height = 6)
cowplot::plot_grid(p1,p2,p3,p4,nrow = 2,align ='hv',labels = c('a','b','c','d'))
dev.off()

#----------------Spatial_pattern ------------------------
directory='./simu_data/diangeju/'

set.seed(123)
#CSR
CSR <- rpoispp(lambda, win = ground, nsim=199)
#Mat
Mat <- rMatClust(kappa = lambda/ave_offspring_per_cluster, scale = radius, 
                 mu = ave_offspring_per_cluster, win = ground, nsim=199)
#Tho
Tho <- rThomas(kappa = lambda/ave_offspring_per_cluster, scale = radius/4,
               mu = ave_offspring_per_cluster, win = ground, nsim=199)
#Hard
Hard <- rHardcore(beta = lambda, R = HC_R, W = ground, nsim=199)
#Str
Str <- rStrauss(beta = lambda, gamma = 0.5, R = HC_R, W = ground, nsim=199)

#choose a format
dat <- coords(rStrauss(beta = lambda, gamma = 0.5, R = HC_R, W = ground, nsim=1))
dat.ppp <- ppp(dat$x, dat$y, window = ground)

#Calculation,CSR
CSRL <- envelope(dat.ppp, Lest, simulate=CSR)
CSRO <- envelope(dat.ppp, pcf, simulate=CSR)

MatL <- envelope(dat.ppp, Lest, simulate=Mat)
MatO <- envelope(dat.ppp, pcf, simulate=Mat)

ThoL <- envelope(dat.ppp, Lest, simulate= Tho)
ThoO <- envelope(dat.ppp, pcf, simulate= Tho)

HCL <- envelope(dat.ppp, Lest, simulate=Hard)
HCO <- envelope(dat.ppp, pcf, simulate=Hard)

StrL<- envelope(dat.ppp, Lest, simulate=Str)
StrO<- envelope(dat.ppp, pcf, simulate=Str)

write.table(CSRL,file = 'simu_data/diangeju/CSR_L.csv',sep = ',',row.names = F)
write.table(CSRO,file = 'simu_data/diangeju/CSR_O.csv',sep = ',',row.names = F)
write.table(MatL,file = 'simu_data/diangeju/Mat_L.csv',sep = ',',row.names = F)
write.table(MatO,file = 'simu_data/diangeju/Mat_O.csv',sep = ',',row.names = F)
write.table(ThoL,file = 'simu_data/diangeju/Tho_L.csv',sep = ',',row.names = F)
write.table(ThoO,file = 'simu_data/diangeju/Tho_O.csv',sep = ',',row.names = F)
write.table(HCL,file = 'simu_data/diangeju/HC_L.csv',sep = ',',row.names = F)
write.table(HCO,file = 'simu_data/diangeju/HC_O.csv',sep = ',',row.names = F)
write.table(StrL,file = 'simu_data/diangeju/Str_L.csv',sep = ',',row.names = F)
write.table(StrO,file = 'simu_data/diangeju/Str_O.csv',sep = ',',row.names = F)

if (F) {
  CSRL<-read.csv('simu_data/diangeju/CSR_L.csv')
  CSRO<-read.csv('simu_data/diangeju/CSR_O.csv')
  MatL<-read.csv('simu_data/diangeju/Mat_L.csv')
  MatO<-read.csv('simu_data/diangeju/Mat_O.csv')
  ThoL<-read.csv('simu_data/diangeju/Tho_L.csv')
  ThoO<-read.csv('simu_data/diangeju/Tho_O.csv')
  HCL<-read.csv('simu_data/diangeju/HC_L.csv')
  HCO<-read.csv('simu_data/diangeju/HC_O.csv')
  StrL<-read.csv('simu_data/diangeju/Str_L.csv')
  StrO<-read.csv('simu_data/diangeju/Str_O.csv')
  data.frame(group=factor(rep(c('Tho','Mat','CSR','Str','HC'),each=nrow(CSRL)),
                          levels = c('Tho','Mat','CSR','Str','HC')),
             rbind(CSRL,MatL,ThoL,HCL,StrL))->all_L
  
  data.frame(group=factor(rep(c('Tho','Mat','CSR','Str','HC'),each=nrow(CSRL)),
                          levels = c('Tho','Mat','CSR','Str','HC')),
             rbind(CSRO,MatO,ThoO,HCO,StrO))->all_O
  
  data.frame(func=rep(c('L','O'),each=nrow(all_L)),rbind(all_L,all_O))->all
  
  write.table(all,file = 'simu_data/diangeju/diangeju_dat.csv',sep = ',',row.names = F)
}
all_L<-read.csv('simu_data/diangeju/diangeju_dat.csv')%>%filter(func=='L')%>%mutate(group=factor(group,levels = c('Tho','Mat','CSR','Str','HC')))
all_O<-read.csv('simu_data/diangeju/diangeju_dat.csv')%>%filter(func=='O')%>%mutate(group=factor(group,levels = c('Tho','Mat','CSR','Str','HC')))

p1<-ggplot(all_L, aes(x=r, y=obs)) +
  geom_ribbon(aes(xmin=0,xmax=50,ymin=lo,ymax=hi,fill=group),alpha=0.5)+
  scale_fill_manual(values =brewer.pal(9,'Set1'))+
  xlab("Scale(m)") + ylab("L(r)") + theme_bw()+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.key.size = unit(15, "pt"),
        legend.position = c(0.1,0.8))

p2<-ggplot(all_O, aes(x=r, y=obs)) +
  geom_ribbon(aes(xmin=min(r),xmax=max(r),ymin=lo,ymax=hi,fill=group),alpha=0.5)+
  scale_fill_manual(values =brewer.pal(9,'Set1'))+
  xlim(c(0,20))+ylim(c(0,15))+
  xlab("Scale(m)") + ylab("g(r)") + theme_bw()+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.key.size = unit(15, "pt"),
        legend.position = c(0.9,0.8))
pdf('figures/fig5.pdf',width = 9,height = 5)
cowplot::plot_grid(p1,p2,align = 'h',labels = letters[1:2])
dev.off()

# ----------------Interpolation --------------------------
directory='./simu_data/interpolation/'
set.seed(1234)
disofd=data.frame();
for (mod in c('Tho','Mat','CSR','Str','HC')){
  cat('=============running ',mod,' ==============\n')
  make_mod(mod = mod)->mod1
  make_net(mod1,type = 'CS')->tree_net
  p<-graph_from_data_frame(tree_net,directed = F)
  
  plants.coords<-coords(mod1)
  node<-unique(c(tree_net$Source,tree_net$Target))
  plants.coords[rownames(plants.coords) %in% node,]->coo
  coo$mdegree<-as.vector(degree(p))
  coo$betweennesssT<-as.vector(betweenness(p,normalized = T))
  coo$betweennesss<-as.vector(betweenness(p))
  coo$closeness<-as.vector(closeness(p))
  write.csv(coo,file = paste0(directory,'/',mod,'_coo.csv'),row.names = F)
}

if(T){
  coo1<-read.csv('simu_data/interpolation/Tho_coo.csv')
  coo2<-read.csv('simu_data/interpolation/Mat_coo.csv')
  coo3<-read.csv('simu_data/interpolation/CSR_coo.csv')
  coo4<-read.csv('simu_data/interpolation/Str_coo.csv')
  coo5<-read.csv('simu_data/interpolation/HC_coo.csv')
}

data.frame(group=factor(rep(c('Tho','Mat','CSR','Str','HC'),
                            c(nrow(coo1),nrow(coo2),nrow(coo3),nrow(coo4),nrow(coo5))),
                        levels = c('Tho','Mat','CSR','Str','HC')),
           rbind(coo1,coo2,coo3,coo4,coo5))->all_coo
write.table(all_coo,file = paste0(directory,'/','interpolation_dat.csv'),sep = ',',row.names = F)
all_coo%>%group_by(group)%>%summarise(mean(mdegree))

if(T){
  Thocha1<- chazhiplot(coo1)
  Matcha1<-  chazhiplot(coo2) 
  CSRcha1<-  chazhiplot(coo3)
  Strcha1<-  chazhiplot(coo4)
  HCcha1<-  chazhiplot(coo5)
  
  Thocha2<- chazhiplot2(coo1) 
  Matcha2<-  chazhiplot2(coo2) 
  CSRcha2<-  chazhiplot2(coo3) 
  Strcha2<-  chazhiplot2(coo4)
  HCcha2<-  chazhiplot2(coo5)
}
pdf(file = 'figures/figS1.pdf',width = 9,height = 16)
#plot_grid(CSRcha1+labs(title = 'CSR'),Matcha1+labs(title = 'Mat'),Thocha1+labs(title = 'Tho'),HCcha1+labs(title = 'HC'),Strcha1+labs(title = 'Str'),rows = 5,align = 'v',labels = 'AUTO')
#plot_grid(CSRcha2+labs(title = 'CSR'),Matcha2+labs(title = 'Mat'),Thocha2+labs(title = 'Tho'),HCcha2+labs(title = 'HC'),Strcha2+labs(title = 'Str'),rows = 5,align = 'v',labels = 'AUTO')
cowplot::plot_grid(Thocha1+labs(title = 'Tho'),Thocha2+labs(title = 'Tho'),
          Matcha1+labs(title = 'Mat'),Matcha2+labs(title = 'Mat'),
          CSRcha1+labs(title = 'CSR'),CSRcha2+labs(title = 'CSR'),
          Strcha1+labs(title = 'Str'),Strcha2+labs(title = 'Str'),
          HCcha1+labs(title = 'HC'), HCcha2+labs(title = 'HC'),
          nrow = 5,ncol = 2,align = 'hv',labels = letters[1:10])

dev.off()

# -------------------Case study -----------------------------
directory='./case_data/'
#1.net analyse
#数据处理,前三列必须x,y,crown_r
dat148<-read.csv('case_data/Lidar_AMS3D_inventory_148_1ha_plots_from_systematic_sampling.csv',skip = 13)
dat148%>%mutate(x=x-min(x),y=y-min(y),crown_r=sqrt(crown_area/pi))%>%select(x,y,crown_r,land_cover_plot_id)->dat148
ggscatter(dat148,'x','y')

dat148%>%filter(land_cover_plot_id=='FS_2')->dat1
summary(dat1)

dat148%>%filter(land_cover_plot_id=='FS_1')->dat3
summary(dat3)

write.csv(dat3,file = 'case_data/example_case.csv',row.names = F)
dat3<-read.csv('case_data/example_case.csv')
to.ppp(dat3)->case_dat
plot_mod(case_dat)

net_analyse(case_dat,res_dir = './case_data/net_analyse_res/',n_simu = 20)

#2.spatial pattern
case_L <- envelope(case_dat, Lest)
case_O <- envelope(case_dat, pcf)
sp1<-ggplot(case_O, aes(x=r, y=obs)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),fill='#4DAF4A',alpha=0.5)+
  geom_line(aes(x=r,y=obs))+
  scale_fill_manual(values =brewer.pal(9,'Set1'))+
  xlab("Scale(m)") + ylab("L(r)") + theme_bw(base_size = 14)+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.key.size = unit(15, "pt"),
        legend.position = c(0.1,0.8))
sp2<-ggplot(case_L, aes(x=r, y=obs)) +
  geom_ribbon(aes(ymin=lo,ymax=hi),fill='#4DAF4A',alpha=0.5)+
  geom_line(aes(x=r,y=obs))+
  scale_fill_manual(values =brewer.pal(9,'Set1'))+
  xlab("Scale(m)") + ylab("g(r)") + theme_bw(base_size = 14)+
  theme(legend.title = element_blank(),
        axis.title.y = element_text(face = "italic"),
        legend.key.size = unit(15, "pt"),
        legend.position = c(0.1,0.8))
cowplot::plot_grid(sp1,sp2,align = 'hv',nrow = 1)%>%ggsave(filename = 'case_data/spa_pattern.pdf',plot = .,width = 7,height = 4)


#3.interpolation
if(T){
CS<-make_net(case_dat,type = 'CS')
p<-graph_from_data_frame(CS,directed = F)
plants.coords<-coords(case_dat)
node<-unique(c(CS$Source,CS$Target))
plants.coords[rownames(plants.coords) %in% node,]->coo
coo$mdegree<-as.vector(degree(p))
coo$betweennesssT<-as.vector(betweenness(p,normalized = T))
coo$betweennesss<-as.vector(betweenness(p))
coo$closeness<-as.vector(closeness(p))
write.csv(coo,file = paste0(directory,'/','CS_case_coo.csv'),row.names = F)
}
#coo<-read.csv('case_data/CS_case_coo.csv')
p1<-chazhiplot(coo)
p2<-chazhiplot2(coo)

cowplot::plot_grid(p1,p2,align = 'hv',nrow = 1)%>%ggsave(filename = 'case_data/interpolation.pdf',plot = .,width = 9,height = 3.2)

nep1<-~{par(mar=c(0,0,1,0),mfcol=c(1,3)); plot_net(case_dat,'CS');plot_net(case_dat,'CL'); plot_net(case_dat,'WCL')}
cowplot::plot_grid(nep1) %>%ggsave(filename = 'rew.pdf',device = 'pdf',width = 9,height = 4)

cowplot::plot_grid(plot_mod(case_dat)[[1]],plot_mod(case_dat)[[2]],sp1,sp2,nrow = 3,labels = 'auto')%>%
  ggsave('figures/fig6.pdf',plot = .,width = 8,height = 10)
#------------Test parameters--------------
#set the parameters
ground <- owin(xrange = c(0, 200), yrange = c(0,200)) #set the investigation area
max_crown <- 5
min_crown <- 2
lambda <- 0.015 #
radius <- 5  #Diffusion radius
ave_offspring_per_cluster <- 3#Mat model
HC_R<-4 #HC model

n_simu=20
directory='./test_pars/'

select_nets=c('CS','CL','WCL')
select_inds=c('mdegree','clusteringC','edensity','aplength')
#-test lambda
set.seed(123)
testpar=data.frame()
for (lambda in c(0.007,0.01,0.013,0.015,0.017,0.02)){
  cat('=============!! test ',lambda,' !!==============\n')
  ttype=data.frame()
  for (type in c('CS','CL','WCL')){
    cat('=============!! running ',type,' !!==============\n')
    tmod=data.frame()
    for (mod in c('Tho','Mat','CSR','Str','HC')){
      cat('=============running ',mod,' ==============\n')
      indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                         mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
      for (j in 1:n_simu){
        make_mod(mod = mod)->mod1
        make_net(mod1,type = type)->tree_net
        #write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
        net_par(tree_net,type = type)->net_indx
        indexs<-rbind(indexs,net_indx)
        if (j%%10==0)cat(j,'have done\n')
      }
      tmod=rbind(tmod,data.frame(mod=mod,num=1:n_simu,indexs))
    }
    ttype=rbind(ttype,data.frame(net=type,tmod))
  }
  testpar=rbind(testpar,data.frame(lambda=lambda,ttype))
}
write.csv(testpar,paste0(directory,'/','test_lambda.csv'),row.names = F)

test_lambda<-read.csv('test_pars/test_lambda.csv')
test_lambda%>%filter(net=='CS',lambda==0.02)->aa
multitest(aa[,'mdegree'],aa[,'mod'])

gglinepls<-list()
for (net in select_nets){
  for (indexs in select_inds){
    test_lambda%>%filter(net==!!net,)%>%mutate(lambda=factor(lambda))%>%rename(indexs=!!indexs)->tmpdf
    tmpp<-ggline(tmpdf,x = 'mod',y ='indexs',color ='lambda',add=c('mean_se'),add.params = list(size=5))+
      ggsci::scale_color_npg()+ylab(label = NULL)+xlab(label = NULL)+theme(legend.position = 'right')
    gglinepls[[paste0(net,'_',indexs)]]<-tmpp
  }
}
pdf('figures/figS2.pdf',14,11)
patchwork::wrap_plots(gglinepls)+  plot_layout(nrow = 3,guides = 'collect')->pppp
pppp+plot_annotation(tag_level = 'a',theme = theme(plot.margin=unit(c(2,0,0,5),'lines')))& 
  theme(plot.tag = element_text(size=22,face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
dev.off()
#-test radius

testpar=data.frame()
for (radius in c(3,4,5,6,7)){
  cat('=============!! test ',radius,' !!==============\n')
  ttype=data.frame()
  for (type in c('CS','CL','WCL')){
    cat('=============!! running ',type,' !!==============\n')
    tmod=data.frame()
    for (mod in c('Tho','Mat','CSR','Str','HC')){
      cat('=============running ',mod,' ==============\n')
      indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                         mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
      for (j in 1:n_simu){
        make_mod(mod = mod)->mod1
        make_net(mod1,type = type)->tree_net
        #write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
        net_par(tree_net,type = type)->net_indx
        indexs<-rbind(indexs,net_indx)
        if (j%%10==0)cat(j,'have done\n')
      }
      tmod=rbind(tmod,data.frame(mod=mod,num=1:n_simu,indexs))
    }
    ttype=rbind(ttype,data.frame(net=type,tmod))
  }
  testpar=rbind(testpar,data.frame(radius=radius,ttype))
}
write.csv(testpar,paste0(directory,'/','test_radius.csv'),row.names = F)
test_radius<-read.csv('test_pars/test_radius.csv')

gglinepls<-list()
for (net in select_nets){
  for (indexs in select_inds){
    test_radius%>%filter(net==!!net,)%>%mutate(radius=factor(radius))%>%rename(indexs=!!indexs)->tmpdf
    tmpp<-ggline(tmpdf,x = 'mod',y ='indexs',color ='radius',add=c('mean_se'),add.params = list(size=5))+
      ggsci::scale_color_npg()+ylab(label = NULL)+xlab(label = NULL)+theme(legend.position = 'right')
    gglinepls[[paste0(net,'_',indexs)]]<-tmpp
  }
}
pdf('figures/figS3.pdf',14,11)
patchwork::wrap_plots(gglinepls)+  plot_layout(nrow = 3,guides = 'collect')->pppp
pppp+plot_annotation(tag_level = 'a',theme = theme(plot.margin=unit(c(2,0,0,5),'lines')))& 
  theme(plot.tag = element_text(size=22,face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
dev.off()
#-test HC_R

testpar=data.frame()
for (HC_R in c(3,4,5,6)){
  cat('=============!! test ',HC_R,' !!==============\n')
  ttype=data.frame()
  for (type in c('CS','CL','WCL')){
    cat('=============!! running ',type,' !!==============\n')
    tmod=data.frame()
    for (mod in c('Tho','Mat','CSR','Str','HC')){
      cat('=============running ',mod,' ==============\n')
      indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                         mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
      for (j in 1:n_simu){
        make_mod(mod = mod)->mod1
        make_net(mod1,type = type)->tree_net
        #write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
        net_par(tree_net,type = type)->net_indx
        indexs<-rbind(indexs,net_indx)
        if (j%%10==0)cat(j,'have done\n')
      }
      tmod=rbind(tmod,data.frame(mod=mod,num=1:n_simu,indexs))
    }
    ttype=rbind(ttype,data.frame(net=type,tmod))
  }
  testpar=rbind(testpar,data.frame(HC_R=HC_R,ttype))
}
write.csv(testpar,paste0(directory,'/','test_HC_R.csv'),row.names = F)
test_HC_R<-read.csv('test_pars/test_HC_R.csv')

gglinepls<-list()
for (net in select_nets){
  for (indexs in select_inds){
    test_HC_R%>%filter(net==!!net,)%>%mutate(HC_R=factor(HC_R))%>%rename(indexs=!!indexs)->tmpdf
    tmpp<-ggline(tmpdf,x = 'mod',y ='indexs',color ='HC_R',add=c('mean_se'),add.params = list(size=5))+
      ggsci::scale_color_npg()+ylab(label = NULL)+xlab(label = NULL)+theme(legend.position = 'right')
    gglinepls[[paste0(net,'_',indexs)]]<-tmpp
  }
}
pdf('figures/figS4.pdf',14,11)
patchwork::wrap_plots(gglinepls)+  plot_layout(nrow = 3,guides = 'collect')->pppp
pppp+plot_annotation(tag_level = 'a',theme = theme(plot.margin=unit(c(2,0,0,5),'lines')))& 
  theme(plot.tag = element_text(size=22,face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
dev.off()

#-test ave_offspring_per_cluster

testpar=data.frame()
for (ave_offspring_per_cluster in c(2,3,4)){
  cat('=============!! test ',ave_offspring_per_cluster,' !!==============\n')
  ttype=data.frame()
  for (type in c('CS','CL','WCL')){
    cat('=============!! running ',type,' !!==============\n')
    tmod=data.frame()
    for (mod in c('Tho','Mat','CSR','Str','HC')){
      cat('=============running ',mod,' ==============\n')
      indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                         mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
      for (j in 1:n_simu){
        make_mod(mod = mod)->mod1
        make_net(mod1,type = type)->tree_net
        #write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
        net_par(tree_net,type = type)->net_indx
        indexs<-rbind(indexs,net_indx)
        if (j%%10==0)cat(j,'have done\n')
      }
      tmod=rbind(tmod,data.frame(mod=mod,num=1:n_simu,indexs))
    }
    ttype=rbind(ttype,data.frame(net=type,tmod))
  }
  testpar=rbind(testpar,data.frame(ave_offspring_per_cluster=ave_offspring_per_cluster,ttype))
}
write.csv(testpar,paste0(directory,'/','test_ave_offspring_per_cluster.csv'),row.names = F)
test_ave_offspring_per_cluster<-read.csv('test_pars/test_ave_offspring_per_cluster.csv')
colnames(test_ave_offspring_per_cluster)[1]<-'offspring'
gglinepls<-list()
for (net in select_nets){
  for (indexs in select_inds){
    test_ave_offspring_per_cluster%>%filter(net==!!net,)%>%mutate(offspring=factor(offspring))%>%rename(indexs=!!indexs)->tmpdf
    tmpp<-ggline(tmpdf,x = 'mod',y ='indexs',color ='offspring',add=c('mean_se'),add.params = list(size=5))+
      ggsci::scale_color_npg()+ylab(label = NULL)+xlab(label = NULL)+theme(legend.position = 'right')
    gglinepls[[paste0(net,'_',indexs)]]<-tmpp
  }
}
pdf('figures/figS5.pdf',14,11)
patchwork::wrap_plots(gglinepls)+  plot_layout(nrow = 3,guides = 'collect')->pppp
pppp+plot_annotation(tag_level = 'a',theme = theme(plot.margin=unit(c(2,0,0,5),'lines')))& 
  theme(plot.tag = element_text(size=22,face = 'bold'),
        legend.title = element_text(size=15),
        legend.text = element_text(size=12))
dev.off()
