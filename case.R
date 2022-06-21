#======== La Selva 生物站9 个 1 公顷森林地块
dat<-read.csv('D:/edgedownload/AMS3D_lidar_points_clouds/Lidar_AMS3D_inventory_9_1ha_plots_over_field_inventory.csv',skip = 11)

dat%>%mutate(x=x-min(x),y=y-min(y),crown_r=sqrt(crown_area/pi))%>%select(x,y,crown_r,idPlot)%>%filter(crown_r>0)->dat
summary(dat)

ggscatter(dat,'x','y')

#先选一个区吧,这个比较好
dat%>%filter(idPlot=='100')->dat1
summary(dat1)
to.ppp(dat1)%>%plot_mod()
#
dat%>%filter(idPlot=='L3')->dat2
summary(dat2)
to.ppp(dat2)%>%plot_mod()

#========这个有好多种类样地
#1.数据处理,前三列必须x,y,crown_r
dat148<-read.csv('case_data/Lidar_AMS3D_inventory_148_1ha_plots_from_systematic_sampling.csv',skip = 13)
dat148%>%mutate(x=x-min(x),y=y-min(y),crown_r=sqrt(crown_area/pi))%>%select(x,y,crown_r,land_cover_plot_id)->dat148
ggscatter(dat148,'x','y')

dat148%>%filter(land_cover_plot_id=='FS_1')->dat3
summary(dat3)

write.csv(dat3,file = 'case_data/example_case.csv',row.names = F)
dat3<-read.csv('case_data/example_case.csv')

net_analyse(dat3,res_dir = './case_data/net_analyse_res2/',n_simu = 5)
