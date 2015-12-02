chromoHBM3D <-function(d, hbm.res, method, args)
{
  merges = hbm.res$merges
  hbm = hbm.res$hbm
  
  ls = length(merges)
  conf = list()
  # iteratively reconstruct starting from the lowest scale up
  for (scale in 1:ls)
  {
    my.scale = merges[[scale]]
    lc = length(my.scale)
    for (cluster in 1:lc)
    {
      indices = my.scale[[cluster]]
      if (length(indices) > 4)
      {
        my.cluster.dist = d[indices,indices]
        
        my.args = list(my.cluster.dist)
        if (length(args) > 0){
          my.args = append(my.args, args, 1) 
        }
        # reconstruct 3d with method of choice
        conf = do.call(method, my.args)
        # update d 
        # but only with the indices that were merged at this scale 
        my.cluster.dist.predict = as.matrix(dist(conf))
        hbm.cluster = hbm[indices,indices]
        dist.to.change = which(hbm.cluster == scale)
        my.cluster.dist[dist.to.change] = my.cluster.dist.predict[dist.to.change]
        d[indices,indices] = my.cluster.dist  
      }
      
    }
  }
  my.args = list(my.cluster.dist)
  if (length(args) > 0){
    my.args = append(my.args, args, 1) 
  }
  conf =  do.call(method, my.args)
  return(conf)
}