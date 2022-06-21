# Forest_net
 Study on the Spatial Structure of Forests from the Perspective of Network  
authors: Peng Chen, Wu Haoran
![](./figures/fig1.png "summary")
## Funtions
1. `net_analyses` can do the whole analyses of a ppp class object come from `make_mod` or case data automatically. All results can be find in a directory, including a map plot, networks plot, all indexs table, and, the `res.pdf`.  
2. `to.ppp` can transfer your coordinate and crown radius data to a ppp class object.   
1. `make_mod` can construct five spatial models('CSR','Mat','HC','Tho','Str').  
2. `plot_mod` can plot the dot map of ppp class object come from `make_mod` or case data.  
3. `make_net` can construct three type of networks ('CS', 'CL', 'WCL') from a ppp class object. its result is a 'source-target-(weight)' dataframe, `graph_from_data_frame` in `igraph` package can handle this.   
4. `plot_net` can plot the three types of networks of ppp class object come from `make_mod` or case data.   
5. `chazhiplot` can do a interpolation on network metrics.

## Main
```
source('funcitons.R')
library(spatstat)
library(igraph)
library(dplyr)
library(reshape2)
library(ggpubr)
library(RColorBrewer)
```

