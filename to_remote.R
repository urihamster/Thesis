require(igraph)
require(Matrix)
require(ggplot2)
require(ROCR)
require(prob)
require(SnowballC)
require(lsa)
require(scales)
library(reshape2)
require(irlba)


datamine_graphs<-function(g1,g2,s)
{
  g2<-g2+vertices(setdiff(V(g1),V(g2)))
  g1<-g1+vertices(setdiff(V(g2),V(g1)))
  
  adj_a<-get.adjacency(g1, type=c("both"),edges=FALSE, names=TRUE,sparse=TRUE)
  adj_b<-get.adjacency(g2, type=c("both"),edges=FALSE, names=TRUE,sparse=TRUE)
  
  adj_a<-adj_a[order(rownames(adj_a)),order(colnames(adj_a))]
  adj_b<-adj_b[order(rownames(adj_b)),order(colnames(adj_b))]
  
  dist_array<-find_distance_vs_new_edges(adj_a,adj_b)
  
  #  svd_a<-irlba(adj_a, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
  #               sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
  #               tol = 1e-06, V = NULL, matmul = NULL)
  
  svd_a<-irlba(adj_a, nu = s, nv = s, maxit = 1000, reorth = 2, tol = 1e-06)
  
  
  #  svd_b<-irlba(adj_b, nu = s, nv = s, adjust = 3, aug = c("ritz","harm"),
  #               sigma = c("ls","ss"), maxit = 1000, m_b = 20, reorth = 2,
  #               tol = 1e-06, V = NULL, matmul = NULL)
  
  svd_b<-irlba(adj_b, nu = s, nv = s, maxit = 1000, reorth = 2, tol = 1e-06)
  
  ne <- adj_b - adj_a
  
  b_svd <- t(svd_a$u) %*% ne %*% svd_a$v
  
  df = as.data.frame(as.numeric(diag(b_svd)))
  
  df$flamda = svd_a$d #the dependent function
  
  colnames(df) = c("source","target")
  
  plot(df)
  
  fit1 <-lm(formula = target ~ I(source) + I(source^2) + I(source^3) + I(source^4) + I(source^5) + I(source^6) + I(source^7) + I(source^8) + I(source^9), data = df)
  
  coefs<-fit1$coefficients
  
  stats_vec1<-get_stats(g1)
  stats_vec2<-get_stats(g2)
  
  the_data<-list(coefs[-1]/sum(coefs[-1]),dist_array)
  
  return(the_data)
}


get_stats<-function(g)
{
  edge_number<-length(E(g))
  m_betweenness<-mean(betweenness(g))
  m_closeness<-mean(closeness(g))
  a_dist<-average.path.length(g)
  stat_vec<-list(edge_number, m_betweenness, m_closeness, a_dist)
  return(stat_vec)
}


find_distance_vs_new_edges<-function(adj1,adj2)
{
  g1<-graph.adjacency(adj1)
  g2<-graph.adjacency(adj2)
  
  dist_ar<-array(1,diameter(g1))
  dist_ar<-dist_ar-dist_ar
  
  dist_mat<-shortest.paths(g1)
  
  diff_mat = adj2-adj1
  
  positives<-which(diff_mat>0,arr.in=TRUE)
  
  for (i in 1:nrow(positives))
  {
    current_dist<-dist_mat[positives[i,][[1]],positives[i,][[2]]]
    dist_ar[current_dist] = dist_ar[current_dist]+1
  }
  
  return(dist_ar)
}



figure_if_edge<-function(g,v1,v2,weights_vec)
{
  dist<-shortest.paths(g,v=v1,to=v2)[1]
  if (dist>1 && !is.na(weights_vec[dist-1]))
  {
    randy<-runif(1)
    if (randy<weights_vec[dist-1])
    {return(1)}
  }
  return(0)
}

generate_next_graph<-function(g,weights_vec,sample1,sample2)
{
  for (i in 1:length(sample1))
  {
    for (j in 1:length(sample2))
    {
      new_edge<-figure_if_edge(g,V(g)[sample1[i]],V(g)[sample2[j]],weights_vec)
      if (new_edge)
      {
        g<-add.edges(g,c(V(g)[sample1[i]],V(g)[sample2[j]]))
      }
    }
  }
  return(g)
}

generate_weights<-function(initial,move)
{
  w_vec<-list(initial) 
  for (i in 2:30)
  {
    w_vec[[i]]<-initial+move*i
  }
  return(w_vec)
}

generate_result<-function(g1,g2,w,s1,s2,s)
{
  sample_vertices1<-sample(V(g2),s1)
  sample_vertices2<-sample(V(g2),s2)
  
  nu_graph<-generate_next_graph(g2,w,sample_vertices1,sample_vertices2)
  
  result<-datamine_graphs(g1,nu_graph,s)
  
  return(result)
}

get_results<-function(g1,g2,the_weights,the_addition,s1,s2,s)
{
  weight_vec<-generate_weights(the_weights,the_addition)
  
  the_results<-list()
  
  for (i in 1:length(weight_vec))
  {
    the_results[[i]]<-generate_result(g1,g2, weight_vec[[i]],s1,s2,s)
  }
  
  return(the_results)
}


g_a<-read.csv('bipartite/uri/as19971108.csv',sep = "\t",header = FALSE, row.names = NULL, col.names = c('as1','as2','relationship'))

edgelist_source <- read.csv(file.choose(),sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))
edgelist_target <- read.csv(file.choose(),sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))
sample_a<-edgelist_source[,1:2] #after data frame
sample_b<-edgelist_target[,1:2] #before data frame
g1<-graph.data.frame(sample_a,directed = FALSE)
g2<-graph.data.frame(sample_b,directed = FALSE)


bla<-1
save.image(file='Thesis/trythisthing.RData')





