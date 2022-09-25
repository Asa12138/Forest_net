#============functions for forest_net analyse============

#Model Construction
make_mod <- function(mod='CSR',crown_r='unif',a=2,b=5,lambda=0.015,ground=owin(xrange = c(0, 200), yrange = c(0,200)),
                       radius=5,HC_R=4,ave_offspring_per_cluster=3) {
  if(!("spatstat" %in% (.packages()))){
    library(spatstat)
  }
  
  if (!mod%in%c('CSR','Mat','HC','Tho','Str'))stop('mod is must be one of CSR,Mat,HC,Tho,Str')
  switch (mod,
          'Tho'={ppp=rThomas(kappa = lambda/ave_offspring_per_cluster, scale = radius/4,
                             mu = ave_offspring_per_cluster, win = ground)},
          'Mat'={ppp=rMatClust(kappa = lambda/ave_offspring_per_cluster, scale = radius, 
                               mu = ave_offspring_per_cluster, win = ground)},
          'CSR'={ppp=rpoispp(lambda, win = ground)},
          'Str'={ppp=rStrauss(beta = lambda, gamma = 0.5, R = HC_R, W = ground)},
          'HC'={ppp=rHardcore(beta = lambda, R = HC_R, W = ground)}
  )
  if (crown_r=='unif')marks(ppp) <- runif(npoints(ppp), max = b, min = a)
  else if (crown_r=='gamma')marks(ppp) <- rgamma(npoints(ppp),shape = a,rate = b)
  
  return(ppp)
}
#Network Construction 
make_net<-function(ppp,type='CS',zb=T,cs_r=10){
  if(!("spatstat" %in% (.packages()))){
    library(spatstat)
  }
  if (!type%in%c('CS','CL','WCL'))stop('type is must be one of CS,CL,WCL')
  plants <- ppp#a ppp class variable
  plants.coords <- coords(plants)#extract the coordinate
  network <- data.frame(Source = c(), Target = c())
  zheban<-\(network){
    na=data.frame()
    as.matrix(network)->a
    for (i in 1:nrow(a)){
      tmp=a[i,]
      na=rbind(na,sort(tmp))
    }
    na[!duplicated(na),]->na
    colnames(na)<-c('Source','Target')
    return(na)
  }
  
  if(type=='CS'){
    for(i in 1:npoints(plants)){
      #Identify overlapped tree
      coords.i <- plants.coords[i,]
      d <- sqrt((plants.coords[,1] - coords.i[,1])^2 + (plants.coords[,2] - coords.i[,2])^2)
      sum_r <- cs_r
      overlap <- which(d < sum_r & d > 0)
      res_i <- data.frame(Source = rep(i,length(overlap)), Target = overlap)
      network <- rbind(network, res_i)
    }
    if(zb)network <-zheban(network)
  }
  if(type=='CL'){
    for(i in 1:npoints(plants)){
      #Identify overlapped tree
      Ri <- marks(plants)[i]
      coords.i <- plants.coords[i,]
      d <- sqrt((plants.coords[,1] - coords.i[,1])^2 + (plants.coords[,2] - coords.i[,2])^2)
      sum_r <- Ri + marks(plants)
      overlap <- which(d < sum_r & d > 0)
      
      res_i <- data.frame(Source = rep(i,length(overlap)), Target = overlap)
      network <- rbind(network, res_i)
    }
    if(zb)network <-zheban(network)
  }
  if(type=='WCL'){
    network <- data.frame(Source = c(), Target = c(),weight =c())
    for(i in 1:npoints(plants)){
      #Identify overlapped tree
      Ri <- marks(plants)[i]
      coords.i <- plants.coords[i,]
      d <- sqrt((plants.coords[,1] - coords.i[,1])^2 + (plants.coords[,2] - coords.i[,2])^2)
      sum_r <- Ri + marks(plants)
      overlap <- which(d < sum_r & d > 0)
      #Calculate HC index
      d.over <- sqrt((plants.coords[overlap,1] - coords.i[,1])^2 +
                       (plants.coords[overlap,2] - coords.i[,2])^2)
      Rj <- marks(plants)[overlap]
      Rij <- (Ri + Rj) - d.over
      CIij <- Rij*Ri/Rj
      
      res_i <- data.frame(Source = rep(i,length(overlap)), Target = overlap,weight = round(CIij,5))
      network <- rbind(network, res_i)
    }
  }
  
  return(network)
}
#Calculate Network Parameters
net_par<-function(net,type='CS'){
  if(!("igraph" %in% (.packages()))){
    library(igraph)
  }
  network<-net
  if (!type%in%c('CS','CL','WCL'))stop('type is must be one of CS,CL,WCL')
  if(type%in%c('CS','CL')){
    p<-graph_from_data_frame(network,directed = F)
    mdegree=mean(degree(p))
  }
  if(type=='WCL'){
    p<-graph_from_data_frame(network)
    mdegree=sum(E(p)$weight)/length(V(p))
  }
  nnodes= length(V(p))#number of nodes
  nedges = length(E(p))#number od edges
  edensity = edge_density(p)#density of network
  aplength=average.path.length(p)
  diameter=diameter(p)
  clusteringC=transitivity(p)
  cb=centralization.betweenness(p)$centralization
  cd=centralization.degree(p)$centralization
  
  index=data.frame(nnodes,nedges,mdegree,clusteringC,edensity,aplength,diameter,cb,cd)
  return(index)
}

plot_net<-function(case_dat,type='CS',size_thr=10,main=NULL){
  if(!("igraph" %in% (.packages()))){
    library(igraph)
  }
  if (!type%in%c('CS','CL','WCL'))stop('type is must be one of CS,CL,WCL')
  make_net(case_dat,type=type)->network
  vertices = data.frame(name=1:case_dat$n,x=case_dat$x,y=case_dat$y,r=case_dat$marks)
  if(is.null(main))main=paste0(type," network")
  if(type%in%c('CS','CL')){
    WCL_net<-graph_from_data_frame(network,vertices = vertices,directed = F)
    V(WCL_net)$size= mmscale(V(WCL_net)$r,5,15)
    V(WCL_net)$color = as.character(cut(degree(WCL_net),breaks =9,labels = brewer.pal(n = 9, name = "Greens")))
    lo <- layout.norm(as.matrix(data.frame(x=case_dat$x,y=case_dat$y)))
    set.seed(123)
    plot(WCL_net,main=main,layout=lo,
         vertex.label.font = 2, vertex.label.color = "black",
         vertex.label.cex = .8,vertex.size = V(WCL_net)$size,
         vertex.label = ifelse(V(WCL_net)$size > size_thr, V(WCL_net)$name, NA),
         edge.lty=1,edge.color='black',edge.curved=F,margin=c(0,0,0,0))
  }
  if(type=='WCL'){
    WCL_net<-graph_from_data_frame(network,vertices = vertices)
    E(WCL_net)$width = mmscale(log2(E(WCL_net)$weight+1),0.5,1)
    V(WCL_net)$size= mmscale(V(WCL_net)$r,5,15)
    V(WCL_net)$color = as.character(cut(degree(WCL_net),breaks =9,labels = brewer.pal(n = 9, name = "Greens")))
    E(WCL_net)$color = as.character(cut(E(WCL_net)$width,breaks =5,labels = brewer.pal(n = 9, name = "YlOrRd")[3:7]))
    
    set.seed(123)
    lo <- layout.norm(as.matrix(data.frame(x=case_dat$x,y=case_dat$y)))
    plot(WCL_net,main=main,layout=lo,
         vertex.label.font = 2, vertex.label.color = "black",
         vertex.label.cex = .8,vertex.size = V(WCL_net)$size,
         vertex.label = ifelse(V(WCL_net)$size > size_thr, V(WCL_net)$name, NA),
         edge.arrow.size=0.2,
         edge.lty=1,edge.curved=F,margin=c(0,0,0,0))
  }
}


net_analyse<-function(case_dat,res_dir='./net_analyse_res/', n_simu=30){
  suppressMessages(lapply(c('ggpubr','igraph','spatstat','fitdistrplus'),library,character.only = TRUE))
  
  #2.casedata构建网络
  res_dir=paste0(res_dir,'/')
  if (!dir.exists(res_dir))dir.create(res_dir)
  #map
  patchwork::wrap_plots(plot_mod(case_dat))%>%ggsave(paste0(res_dir,'case_map.pdf'),plot = .,width = 12,height = 6)
  cat('1.data input& map plot done \n')
  #case_net
  WCL<-make_net(case_dat,type = 'WCL')
  CL<-make_net(case_dat,type = 'CL')
  CS<-make_net(case_dat,type = 'CS')
  netindexs<-rbind(net_par(CS,type = 'CS'),net_par(CL,type = 'CL'),net_par(WCL,type = 'WCL'))
  rownames(netindexs)<-c('CS','CL','WCL')
  write.csv(netindexs,file = paste0(res_dir,'case_netindexs.csv'))
  
  write.table(WCL,file = paste0(res_dir,'case_WCL.csv'),row.names = F,sep = ',')
  write.table(CL,file = paste0(res_dir,'case_CL.csv'),row.names = F,sep = ',')
  write.table(CS,file = paste0(res_dir,'case_CS.csv'),row.names = F,sep = ',')
  #plot network
  pdf(paste0(res_dir,'case_networks.pdf'),width = 12,height = 5)
  par(mfcol=c(1,3)) 
  plot_net(case_dat,'CS')
  plot_net(case_dat,'CL')
  plot_net(case_dat,'WCL')
  dev.off()
  cat('2.network construct& plot done \n')
  #3.构建null model且比对
  #这个样地的case树冠比较符合一个gamma分布

  fit <- fitdist(case_dat$marks, distr = "gamma", method = "mle")
  shape=fit$estimate[1];rate=fit$estimate[2]
  fit_marks=rgamma(length(case_dat$marks),shape,rate)
  rbind(data.frame(type='case',marks=case_dat$marks),
        data.frame(type='fit',marks=fit_marks))%>%
    ggdensity(x = 'marks',color='type',title = 'Crown_R')%>%
    ggsave(paste0(res_dir,'/','crown_R_fit.pdf'),plot = .,width = 6,height = 5)
  
  #set the parameters
  ground <- case_dat$window
  lambda <-  case_dat$n/(ground$xrange*ground$yrange)[2]
  radius <- 5  #Diffusion radius
  HC_R<-4 #HC model
  ave_offspring_per_cluster <- 3
  cat('3.set parameters done& simulation start \n')
  #设置参与并行的CPU核数目，这里我们使用了所有的CPU核，也就是我们刚才得到的clnum，具体到这个案例，clnum=8
  no_cores<-detectCores(logical=F)
  cl <- makeCluster(no_cores);
  #simulate
  set.seed(1234)
  pb=tpb(0,3*5,style = 3);cnt=1
  
  ttype=data.frame()
  for (type in c('CS','CL','WCL')){
    #cat('=============!! running ',type,' !!==============\n')
    #dir.create(paste0(directory,type))
    tmod=data.frame()
    for (mod in c('Tho','Mat','CSR','Str','HC')){
      #cat('=============running ',mod,' ==============\n')
      indexs<-data.frame(nnodes= c(), nedges = c(), edensity = c(),aplength=c(),
                         mdegree=c(),diameter=c(),clusteringC=c(),cb=c(),cd=c())
      processbar(pb,cnt)
      fun<-\(x){
        make_mod(mod = mod,crown_r = 'gamma',a = shape,b=rate,ground = ground,lambda = lambda,
                 radius = radius,HC_R = HC_R,ave_offspring_per_cluster = ave_offspring_per_cluster)->mod1
        make_net(mod1,type = type)->tree_net
        #write.csv(tree_net,file = paste0(paste0(directory,type),'/',mod,'_',j,'.csv'),row.names = F)
        net_par(tree_net,type = type)->net_indx
        #if (j%%10==0)cat(j,'have done\n')
        (net_indx%>%t())[,1]->net_indx
        return(net_indx) 
      }
      clusterExport(cl, c('make_mod','make_net','net_par','mod','type','shape','rate','ground','lambda','radius','HC_R','ave_offspring_per_cluster'))
      indexs<-parSapply(cl,1:n_simu,fun)%>%t()
      
      tmod=rbind(tmod,data.frame(mod=mod,num=1:n_simu,indexs))
      cnt=cnt+1
    }
    ttype=rbind(ttype,data.frame(net=type,tmod))
  }
  #关闭并行计算
  stopCluster(cl);
  
  write.csv(ttype,file = paste0(res_dir,'/','null_model_index.csv'),row.names = F)
  
  #statistic
  options(dplyr.summarise.inform = FALSE)
  all_index<-ttype
  all_index%>%group_by(net,mod)%>%summarise(nnodes_min=round(min(nnodes),0), nnodes_mean=round(mean(nnodes),0),nnodes_max=round(max(nnodes),0),
                                            nedges_min=round(min(nedges),0),nedges_mean=round(mean(nedges),0),nedges_max=round(max(nedges),0),
                                            k_min=round(min(mdegree),2),k_mean=round(mean(mdegree),2),k_max=round(max(mdegree),2),
                                            C_min=round(min(clusteringC),2),C_mean=round(mean(clusteringC),2),C_max=round(max(clusteringC),2),
                                            D_min=round(min(edensity),4),D_mean=round(mean(edensity),4),D_max=round(max(edensity),4),
                                            L_min=round(min(aplength),2),L_mean=round(mean(aplength),2),L_max=round(max(aplength),2))%>%
    write.csv(file = paste0(res_dir,'/','null_model_stat.csv'),row.names = F)
  
  cat('\n4.null model simulation done \n')
  all_index$net<-factor(all_index$net,levels = c('CS','CL','WCL'))
  all_index$mod<-factor(all_index$mod,levels = c('Tho','Mat','CSR','Str','HC'))
  #all indexs in five models of three networks
  melt(all_index,id.vars = c('net','mod','num'),variable.name = 'indexs')->long_index
  #4.最终结果
  select_nets=c('CS','CL','WCL')
  select_inds=c('mdegree','clusteringC','edensity','aplength')
  boxpls<-list()
  for (net in select_nets){
    for (index in select_inds){
      index1<-switch(index,
                     'mdegree'='k',
                     'clusteringC'='C',
                     'edensity'='D',
                     'aplength'='L')
      boxpls[[paste0(net,'_',index)]]<-box_plot(long_index,net,index)+geom_hline(yintercept = netindexs[net,index],col="#10C0E3",lty=2)+
        annotate(geom = 'text',x = 'Str',y=netindexs[net,index]*1.05,label='case',col="#10C0E3")+labs(title =paste0(net,'_',index1))
    }
  }
  cowplot::plot_grid(plotlist = boxpls,nrow = 3,align ='hv',labels ='auto',label_size = 21)->res_p
  pdf(paste0(res_dir,'/','res.pdf'),16,12)
  plot(res_p)
  dev.off()
  cat('5.result plot done \n')
  cat('Thanks \n')
}

plot_mod<-function(mod1){
  library(ggpubr)
  pld<-data.frame(x=mod1$x,y=mod1$y,crown=mod1$marks)
  p1<-ggscatter(pld,'x','y',size = 'crown',color='green4')+
    geom_density2d(size = 0.5)+labs(x='X/m',y='Y/m')
  
  p2<-ggplot(pld,aes(x,y,size = crown)) +
    ggpointdensity::geom_pointdensity(adjust = 100) +
    viridis::scale_color_viridis()+theme_pubr()+labs(x='X/m',y='Y/m')+guides(size='none')
  return(list(p1,p2))
}

multitest<-\(var,group){
  library(agricolae)
  group<-factor(group)
  aov(var~group)->ano
  #print(summary(ano))
  lsdres <- LSD.test(ano, 'group', p.adj = 'bonferroni')
  data.frame(lsdres$groups)->a
  row.names(a)->a$variable
  #print(TukeyHSD(ano, 'group', p.adj = 'bonferroni'))
  return(a)
}

mmscale=\(x,min=0,max=1){
  return(min+(x-min(x))/(max(x)-min(x))*(max-min))
}

box_plot<-function(long_index,net,indexs){
  if(!("ggpubr" %in% (.packages()))){
    library(ggpubr)
  }
  long_index%>%filter(net=={!!net},indexs=={!!indexs})->tmpdf
  1.1*max(tmpdf$value)->high;
  multitest(tmpdf$value,tmpdf$mod)->lsd1
  tmpp<-ggboxplot(tmpdf,x = 'mod',y='value',color='mod',add = c('jitter'),
                  add.params = list(color='mod',width = 0.15,alpha=0.8,size=0.5),outlier.shape = NA)+
    scale_color_manual(values =brewer.pal(9,'Set1'))+ylab(label = NULL)+xlab(label = NULL)+
    geom_text(data = lsd1,aes(x=variable,y=high,label=groups), inherit.aes = FALSE ,color='red',size=5)+
    theme_pubr(base_size = 14)+ theme(plot.margin=unit(rep(1,4),'lines'),legend.position="none")
  tmpp
}

displot<-function(dis){
  library(ggridges)
  library(viridis)
  library(hrbrthemes)
  colnames(dis)<-c('mod','degree')
  dis$mod<-factor(dis$mod,levels = rev(c('Tho','Mat','CSR','Str','HC')))
  p<-ggplot(dis, aes(x =degree , y = mod, fill = ..x..)) +
    geom_density_ridges_gradient(scale = 2, rel_min_height = 0.001) +
    scale_fill_viridis(name = "Temp. [F]", option = "C") +ylab(label = '')+
    theme(
      legend.position="none",
      panel.spacing = unit(0.1, "lines"),
      axis.title.x = element_text(size=12)
    ) + theme(panel.grid.major = element_line(colour = "gray75",size = 0.4),
              panel.grid.minor = element_line(colour = "gray75"),
              panel.background = element_rect(fill = NA),
              plot.title = element_text(size = 18,
                                        face = "bold")) +labs(y = NULL)
  return(p)
}

to.ppp<-function(dat){
  library(spatstat)
  if (!all(tolower(colnames(dat)[1:3])==c('x','y','crown_r')))stop('x,y,crown_r must in the first three column')
  xrange=c(min(dat$x)-min(dat$x)%%5,max(dat$x)+5-max(dat$x)%%5)
  yrange=c(min(dat$y)-min(dat$y)%%5,max(dat$y)+5-max(dat$y)%%5)
  ppp(dat$x-xrange[1],dat$y-yrange[1],window =owin(xrange = xrange-xrange[1], yrange = yrange-yrange[1]), marks = dat$crown_r)->ppp
  return(ppp)
}

toXY <- function(lat, long){
  if(!("SoDA" %in% (.packages()))){
    library(SoDA)
  }
  XY <- geoXY(lat,long)
  return(as.data.frame(XY))
}

geo.interpolate <- function(dat, var, method = c("Spline", "Inverse"), 
                            Xextend = "defaulted", Yextend = "defaulted"){
  # 'spline' 样条函数插值法
  # 'inverse' 逆距离加权插值法
  if(!("sp" %in% (.packages()))){
    library(sp)
  }
  
  if(!("gstat" %in% (.packages()))){
    library(gstat)
  }
  
  if(!("raster" %in% (.packages()))){
    library(raster)
  }
  
  if(!("fields" %in% (.packages()))){
    library(fields)
  }
  
  if(class(dat) != "data.frame"){
    stop("`dat` should be a 'data.frame'")
  }
  
  if(!("X" %in% colnames(dat))){
    stop("can't find 'X' field that specifies longtitude in `dat`. 
         Have you transformed geographic information? Use `toXY` function")
  }
  
  if(!("Y" %in% colnames(dat))){
    stop("can't find 'X' field that specifies longtitude in `dat`. 
         Have you transformed geographic information? Use `toXY` function")
  }
  
  Xmax <- max(dat$X)
  Ymax <- max(dat$Y)
  Xmin <- min(dat$X)
  Ymin <- min(dat$Y)
  if(Xextend == "defaulted"){
    Xextend = (Xmax - Xmin) * 0.1
  }
  if(Yextend == "defaulted"){
    Yextend = (Ymax - Ymin) * 0.1
  }
  
  pol <- Polygon(data.frame(c(Xmin - Xextend,Xmin - Xextend,Xmax + Xextend,Xmax + Xextend,Xmin - Xextend),
                            c(Ymin - Yextend,Ymax + Yextend,Ymax + Yextend,Ymin - Yextend,Ymin - Yextend)))
  Sr <- SpatialPolygons(list(Polygons(list(pol),1)))
  bound <- SpatialPolygonsDataFrame(Sr, data = data.frame(Id = 0))
  
  Co_inf <- SpatialPoints(dat[,1:2]) 
  Co_inf <- SpatialPointsDataFrame(Co_inf,dat)
  blank_raster <- raster(nrow=100,ncol=100,extent(bound)) 
  values(blank_raster) <- 1 
  bound_raster<-rasterize(bound,blank_raster)
  
  if(method == "Spline"){
    m <- Tps(coordinates(Co_inf), var)
  } else if(method == "Inverse"){
    m <- gstat(formula=var~1, locations=Co_inf) 
  } 
  res <- interpolate(bound_raster, m)
  return(res)
}

krige.interpolate <- function(dat, var, psill=NA, model="Sph", nugget=NA, nmax = Inf, nmin = 0,
                              Xextend = "defaulted", Yextend = "defaulted"){
  if(!("sp" %in% (.packages()))){
    library(sp)
  }
  
  if(!("gstat" %in% (.packages()))){
    library(gstat)
  }
  
  if(!("raster" %in% (.packages()))){
    library(raster)
  }
  
  if(!("fields" %in% (.packages()))){
    library(fields)
  }
  
  if(class(dat) != "data.frame"){
    stop("`dat` should be a 'data.frame'")
  }
  
  if(!("X" %in% colnames(dat))){
    stop("can't find 'X' field that specifies longtitude in `dat`. 
         Have you transformed geographic information? Use `toXY` function")
  }
  
  if(!("Y" %in% colnames(dat))){
    stop("can't find 'X' field that specifies longtitude in `dat`. 
         Have you transformed geographic information? Use `toXY` function")
  }
  
  Xmax <- max(dat$X)
  Ymax <- max(dat$Y)
  Xmin <- min(dat$X)
  Ymin <- min(dat$Y)
  if(Xextend == "defaulted"){
    Xextend = (Xmax - Xmin) * 0.1
  }
  if(Yextend == "defaulted"){
    Yextend = (Ymax - Ymin) * 0.1
  }
  
  pol <- Polygon(data.frame(c(Xmin - Xextend,Xmin - Xextend,Xmax + Xextend,Xmax + Xextend,Xmin - Xextend),
                            c(Ymin - Yextend,Ymax + Yextend,Ymax + Yextend,Ymin - Yextend,Ymin - Yextend)))
  Sr <- SpatialPolygons(list(Polygons(list(pol),1)))
  bound <- SpatialPolygonsDataFrame(Sr, data = data.frame(Id = 0))
  
  Co_inf <- SpatialPoints(dat[,1:2]) 
  Co_inf <- SpatialPointsDataFrame(Co_inf,dat)
  blank_raster <- raster(nrow=100,ncol=100,extent(bound)) 
  values(blank_raster) <- 1 
  bound_raster<-rasterize(bound,blank_raster)
  
  if(any(var<0)){
    stop("negative values found in `var`")
  }
  v <- variogram(log(var) ~ 1, data = Co_inf) 
  v.fit<-fit.variogram(v,model=vgm(psill, model, nugget))
  
  Grid <- as(bound_raster,"SpatialGridDataFrame") 
  res <- krige(formula= var~1 ,model=v.fit,locations=Co_inf,newdata=Grid, nmax=nmax, nmin=nmin)   
  
  return(list(variogram = v.fit, krigging = res))
}

as.data.frame.RasterLayer <- function(ras){
  xstep <- (xmax(ras) - xmin(ras))/(nrow(ras)-1)
  ystep <- (ymax(ras) - ymin(ras))/(ncol(ras)-1)
  x <- seq(xmin(ras), xmax(ras), xstep)
  y <- seq(ymin(ras), ymax(ras), ystep)
  
  y <- rev(y) # Transformation
  y.grid <- c()
  for(i in 1:nrow(ras)){ y.grid <- c(y.grid,rep(y[i],nrow(ras))) }
  x.grid <- rep(x,ncol(ras))
  return(data.frame(x=x.grid,y=y.grid,value=values(ras)))
}

chazhiplot<-function(dat){
  library(ggplot2)
  colnames(dat)[1:2]<-c('X','Y')
  ras <- geo.interpolate(dat,dat$mdegree,method = "Spline") #modi
  to.plot <- as.data.frame.RasterLayer(ras)
  g<-ggplot(to.plot, aes(x, y))+ geom_tile(aes(fill = value))  +  xlab("X/m") + ylab("Y/m")+theme(
    legend.key.size = unit(10, "pt") )+ scale_fill_gradient(name="mdegree", low = "#3EF809", high = "#F15716",limits=c(-5,25))
  return(g)
}

chazhiplot2<-function(dat){
  library(ggplot2)
  colnames(dat)[1:2]<-c('X','Y')
  ras <- geo.interpolate(dat,dat$betweennesssT,method = "Spline") #modi
  to.plot <- as.data.frame.RasterLayer(ras)
  g<-ggplot(to.plot, aes(x, y))+ geom_tile(aes(fill = value))  +  xlab("X/m") + ylab("Y/m")+theme(
    legend.key.size = unit(10, "pt") ) + scale_fill_gradient(name="betweennesss", low = "#3EF809", high = "#F15716")
  return(g)
}

processbar<-function (pb, value, title = NULL, label = NULL) 
{
  if (!inherits(pb, "txtProgressBar")) 
    stop(gettextf("'pb' is not from class %s", dQuote("txtProgressBar")), 
         domain = NA)
  oldval <- pb$getVal()
  pb$up(value)
  invisible(oldval)
}
tpb<-function (min = 0, max = 1, initial = 0, char = "=", width = NA, 
               title, label, style = 1, file = "") 
{
  if (!identical(file, "") && !(inherits(file, "connection") && 
                                isOpen(file))) 
    stop("'file' must be \"\" or an open connection object")
  if (!style %in% 1L:3L) 
    style <- 1
  .val <- initial
  .killed <- FALSE
  .nb <- 0L
  .pc <- -1L
  nw <- nchar(char, "w")
  if (nw == 0) 
    stop("'char' must have a non-zero width")
  if (is.na(width)) {
    width <- getOption("width")
    if (style == 3L) 
      width <- width - 10L
    if (nw > 1) 
      width <- trunc(width/nw)
  }
  if (max <= min) 
    stop("must have 'max' > 'min'")
  up1 <- \(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb < nb) {
      cat(strrep(char, nb - .nb), file = file)
      flush.console()
    }
    else if (.nb > nb) {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up2 <- \(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    if (.nb <= nb) {
      cat("\r", strrep(char, nb), sep = "", file = file)
      flush.console()
    }
    else {
      cat("\r", strrep(" ", .nb * nw), "\r", strrep(char, 
                                                    nb), sep = "", file = file)
      flush.console()
    }
    .nb <<- nb
  }
  up3 <- \(value) {
    if (!is.finite(value) || value < min || value > max) 
      return()
    .val <<- value
    nb <- round(width * (value - min)/(max - min))
    pc <- round(100 * (value - min)/(max - min))
    if (nb == .nb && pc == .pc) 
      return()
    cat(paste0("\r  |", strrep(" ", nw * width + 6)), file = file)
    cat(paste(c("\r  |", rep.int(char, nb), rep.int(" ", 
                                                    nw * (width - nb)), sprintf("| %3d%%", pc)), collapse = ""), 
        file = file)
    flush.console()
    .nb <<- nb
    .pc <<- pc
  }
  getVal <- function() .val
  kill <- function() if (!.killed) {
    cat("\n", file = file)
    flush.console()
    .killed <<- TRUE
  }
  up <- switch(style, up1, up2, up3)
  up(initial)
  structure(list(getVal = getVal, up = up, kill = kill), class = "txtProgressBar")
}

