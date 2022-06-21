#============functions for forest_net analyses============

#Model Construction
make_mod <- function(mod='CSR',crown_r='unif') {
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
  if (crown_r=='unif')marks(ppp) <- runif(npoints(ppp), max = max_crown, min = min_crown)
  else if (crown_r=='gamma')marks(ppp) <- rgamma(npoints(ppp),shape = shape,rate = rate)
  
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
      CIij <- Rij*Rj/Ri
      
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

plot_mod<-function(mod1){
  library(ggpubr)
  pld<-data.frame(x=mod1$x,y=mod1$y,crown=mod1$marks)
  p1<-ggscatter(pld,'x','y',size = 'crown',color='green4')+
    geom_density2d(size = 0.5)+labs(x='X/m',y='Y/m')
  
  p2<-ggplot(pld,aes(x,y,size = crown)) +
    ggpointdensity::geom_pointdensity(adjust = 300) +
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
  g<-ggplot(to.plot, aes(x, y))+ geom_tile(aes(fill = value))  +  xlab("X/m") + ylab("Y/m")+ theme_pubr()+theme(
    legend.key.size = unit(10, "pt") )+ scale_fill_gradient(name="mdegree", low = "#3EF809", high = "#F15716")
  return(g)
}

chazhiplot2<-function(dat){
  library(ggplot2)
  colnames(dat)[1:2]<-c('X','Y')
  ras <- geo.interpolate(dat,dat$betweennesssT,method = "Spline") #modi
  to.plot <- as.data.frame.RasterLayer(ras)
  g<-ggplot(to.plot, aes(x, y))+ geom_tile(aes(fill = value))  +  xlab("X/m") + ylab("Y/m")+ theme_pubr()+theme(
    legend.key.size = unit(10, "pt") ) + scale_fill_gradient(name="betweennesss", low = "#3EF809", high = "#F15716")
  return(g)
}