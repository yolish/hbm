library(igraph) # make sure you have version installed >= 1.0.1

chromoHBM <- function(A,  cluster.method, cluster.args)
{
  if (!is.matrix(A))
  {
    stop("first argument must be of type matrix")
  }
  N = nrow(A)
  if (N != ncol(A))
  {
    stop("input matrix must be a square matrix")
  }
  hbm = matrix(0, N, N)
  scales = list()
  scale = 1
  Ncurr = N
  # remove self interactions
  diag(A) = 0
  
  # init the local contact matrix
  cp = A
  
  # iteratively cluster
  while ( ( ((N-Ncurr) > 0) || (length(scales) == 0) ) && (Ncurr > 1))
  {
    N = Ncurr    
    g = graph.adjacency(cp, mode = "undirected", weighted = TRUE)
    
    # cluster with clustering method of choice
    # put the graph to be the first argument
    my.args = list(g)
    if (length(cluster.args) > 0){
      my.args = append(my.args, cluster.args, 1) 
      
    }
    clusters = do.call(cluster.method, my.args)
    
    # clusters should be returned in 1 of 2 forms
    # communities object with membership or as a vector 
    if (class(clusters) == "communities")
    {
      clusters = as.vector(membership(clusters)) # make it a vector
    }
    
    clustid = unique(clusters)
    ids = lapply(clustid, function(id) which(clusters == id))
    Ncurr = length(clustid)
    if (N == Ncurr)
    {
      # cannot merge any more
      ids = list(c(1:N))
    }
    
    # update scales
    scales[[scale]] = ids
    
    if (scale > 1)
    {
      currscale = scales[[scale]]
      prevscale = scales[[(scale-1)]]  
      for (i in 1:length(currscale))
      {
        indices = currscale[[i]]
        orig.ids = unlist(lapply(indices, function(x) prevscale[[x]]))
        currscale[[i]] = orig.ids
      }
      
      scales[[scale]] = currscale
    }
    
    # update the hbm and cp matrices 
    if (!(Ncurr == 1 || length(scales[[scale]]) == 1)) # we're done
    {
      cp = matrix(0, Ncurr, Ncurr)
      for (i in 1:(Ncurr))
      {
        indices.i = scales[[scale]][[i]]
        hbm[indices.i, indices.i] = hbm[indices.i, indices.i] + 1
        for (j in i:Ncurr)
        {
          indices.j = scales[[scale]][[j]]
          cp[i,j] = sum(A[indices.i, indices.j])/(length(indices.j)*length(indices.i))
        }
      }
      diag(cp) = 0
      cp[lower.tri(cp)] = t(cp)[lower.tri(t(cp))] 
    }
    scale = scale+1 
  }
  # reverse scales for visualization and easier comprehension 
  # nodes that were clustered at the ith iteration will correspond to scale entry i 
  max.scale = max(hbm) + 1
  hbm = hbm*(-1) + max.scale
  diag(hbm) = 1
  return(list("hbm" = hbm, "merges" = scales)) 
}







