
# h1 is a mxm hbm and h2 is a kxk hbm, 
# and m = nk where n is a natural number

chromoHBM3C <-function(h1, h2) 
{
    m = nrow(h1)
    k = nrow(h2)
    n = k/m
    
    # expand H1
    h1.expand = matrix(0, k, k)
    diag(h1.expand) = 1
    for (i in 1:(k-1))
    {
      p = ceiling(i/n)
      for (j in (i+1):(k))
      {
        
        q = ceiling(j/n)
        h1.expand[i,j]= h1[p,q]
        h1.expand[j,i] = h1.expand[i,j]
      }
    }
    
    # compute the merging matrix
    merging = (h1.expand + h2)/2
    h.merge = matrix(0, k, k)
    my.levels = sort(unique(as.vector(merging)))
    nlevels = length(my.levels)
    diag(h.merge) = 1
    for (l in 1:nlevels)
    {
      h.merge[which(merging == my.levels[l])] = l
    }
    
    return(h.merge)
}






