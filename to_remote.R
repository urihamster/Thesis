install.packages("igraph")
install.packages("Matrix")
install.packages("ggplot2")
install.packages("ROCR")
install.packages("prob")
install.packages("SnowballC")
install.packages("lsa")
install.packages("scales")
install.packages("reshape2")
install.packages("irlba")


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
    print(i)
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

get_results<-function(g1,g2,the_weights,the_addition,s1,s2,s,the_start, the_end )
{
  weight_vec<-generate_weights(the_weights,the_addition)
  
  the_results<-list()
  
  for (i in the_start:the_end)
  {
    print('new weight!!!!! ---- !!!!! --- !!!')
    print(i)
    print('new weight!!!!! ---- !!!!! --- !!!')
    print('also - the start:')
    print(the_start)
    print('and the end:')
    print(the_end)
    print('hihihihihihihihihihihihihihihihihihihi')
    the_results[[i]]<-generate_result(g1,g2, weight_vec[[i]],s1,s2,s)
  }
  
  return(the_results)
}


mine_g<-function(g1,g2,the_weights,the_addition,s1,s2,s,file_way,the_start, the_end)
{
  results<-get_results(g1,g2,the_weights,the_addition,s1,s2,s, the_start, the_end)
  save(results,file=file_way)
  return(results)
}


edgelist_source <- read.csv('Thesis/as19990427.csv',sep = ",",header=FALSE, row.names = NULL, col.names = c('from','to', 'relation'))
edgelist_target <- read.csv('Thesis/as19990428.csv',sep = ",",header=FALSE, row.names = NULL, col.names = c('from','to', 'relation'))

sample_a<-edgelist_source[,1:2] #after data frame
sample_b<-edgelist_target[,1:2] #before data frame
g1<-graph.data.frame(sample_a,directed = FALSE)
g2<-graph.data.frame(sample_b,directed = FALSE)

w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
a<-c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)

# - you are at 13 and 14 goin on on on on 

results_1<-mine_g(g1,g2,w,a,100,50,40,'Thesis/resultsB24_1.rda',1,10)

results_2<-mine_g(g1,g2,w,a,100,50,40,'Thesis/resultsB24_2.rda',11,20)

results_3<-mine_g(g1,g2,w,a,100,50,40,'Thesis/resultsB24_3.rda',21,30)


edgelist_source <- read.csv('Thesis/as-caida1.csv',sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))

edgelist_source <- read.csv('Thesis/as19971108.csv',sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))
edgelist_target<- read.csv('Thesis/as19971109.csv',sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))


edgelist_source <- read.csv(file.choose(),sep = ",",header=FALSE, row.names = NULL, col.names = c('from','to'))
hmmm <-read.csv(file.choose(),sep = ",",header=FALSE, row.names = NULL, col.names = c('from','to'),skip=3)




bla<-1
save.image(file='Thesis/trythisthing.RData')

g_a<-read.csv('bipartite/uri/as19971108.csv',sep = "\t",header = FALSE, row.names = NULL, col.names = c('as1','as2','relationship'))

edgelist_source <- read.csv(file.choose(),sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))
edgelist_target <- read.csv(file.choose(),sep = ",",header=TRUE, row.names = NULL, col.names = c('from','to'))


results3<-results

for (i in 11:20)
{
  results3[[i]]<-results2[[i]]
}

resultsW1<-results3

results3




# w<-c(0.1, 0.05)
# a<-c(0.02, 0.005)
# B --- Done
resultsW1 # Done!

# w<-c(0.1, 0.05, 0.01)
# a<-c(0.02, 0.005, 0.001)
# results2_1, results2_2  and results2_3
# B --- Done
resultsW2 # Not Done - file done

# w<-c(0.1, 0.05, 0.01)
# a<-c(0.04, 0.01, 0.002)
# results3_1 ...
# B --- DONE
resultsW3 # not Done - file done

#w<-c(0.1, 0.05, 0.01)
#a<-c(0.08, 0.02, 0.004)
# results4_1 ...
# B --- DONE
resultsW4 # not done - file done


#w<-c(0.1, 0.05, 0.01, 0.01)
#a<-c(0.04, 0.01, 0.002, 0.002)
# results5_1 ...
# B ---  DONE
resultsW5 # not done - file done


#w<-c(0.1, 0.05, 0.01, 0.01, 0.01)
#a<-c(0.04, 0.01, 0.002, 0.002, 0.002)
# results6_1 ...
# B ---  DONE
resultsW6 # not done - file done

#w<-c(0.1, 0.05, 0.02, 0.02, 0.02)
#a<-c(0.04, 0.01, 0.004, 0.004, 0.004)
# results7_1 ...
# B --- DONE
resultsW7 # not done - file done

#w<-c(0.1, 0.05, 0.05, 0.05, 0.05)
#a<-c(0.04, 0.01, 0.01, 0.01, 0.01)
# results8_1 ...
# B --- DONE
resultsW8 # not done - file done


#w<-c(0.1, 0.075, 0.075, 0.075, 0.075)
#a<-c(0.04, 0.025, 0.025, 0.025, 0.025)
# results9_1 ...
# B ---  DONE
resultsW9 # not done - file done


#w<-c(0.1, 0.1, 0.1, 0.1, 0.1)
#a<-c(0.04, 0.04, 0.04, 0.04, 0.04)
# results10_1 ...
# B --- DONE
resultsW10 # not done - file done


#w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
#a<-c(0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04, 0.04)
# results11_1 ...
# B ---  DONE
resultsW11 # not done - file done

#w<-c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5)
#a<-c(0, 0, 0, 0, 0, 0, 0, 0, 0)
# results12_1 ...
# B --- DONE
resultsW12 # not done - file done


# w<-c(0.1)
# a<-c(0.03)
# results13_1 ...
# B --- DONE
resultsW13 # not done - file done

# w<-c(0,0.1)
# a<-c(0,0.03)
# results14_1 ...
# B --- DONE
resultsW14 # not done - _1 done. _2 and _3 not yet done


# w<-c(0,0,0.1)
# a<-c(0,0,0.03)
# results15_1 ...
# B --- DONE
resultsW15 # not done 

# w<-c(0,0,0,0.1)
# a<-c(0,0,0,0.03)
# results16_1 ...
# B --- DONE
resultsW16 # not done 

# w<-c(0,0,0,0,0.1)
# a<-c(0,0,0,0,0.03)
# results17_1 ...
# B --- DONE
resultsW17 # not done 

# w<-c(0.1, 0.1)
# a<-c(0.03, 0.03)
# results18_1 ...
# B --- DONE
resultsW18 # not done

# w<-c(0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03)
# results19_1 ...
# B --- DONE
resultsW19 # not done - file done

# w<-c(0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03)
# results20_1 ...
# B --- DONE
resultsW20 # not done - file done

# w<-c(0.1, 0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03, 0.03)
# results21_1 ...
# B --- DONE
resultsW21 # not done - file done

#w<-c(1, 1, 1, 1, 1, 1, 1, 1, 1)
#a<-c(0, 0, 0, 0, 0, 0, 0, 0, 0)
# results22_1 ...
# B --- DONE
resultsW22 # not done - file done

# w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
# results23_1 ...
# B --- DONE
resultsW23 # not done 

# w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
# results24_1 ...
# B --- DONE
resultsW24 # not done 

# w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
# results25_1 ...
# B --- NOT DONE
resultsW25 # not done 


# w<-c(0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1)
# a<-c(0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03, 0.03)
# results26_1 ...
# B --- NOT DONE
resultsW26 # not done 
