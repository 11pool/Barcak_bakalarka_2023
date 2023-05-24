#Kód v jazyku R použitý pre bakalársku prácu:
#Urèovanie poètu zhlukov v dátach
#Adam Barèák
#Univerzita Komenského FMFI, 2023
#####################################
library("cluster")
library("dbscan")
library("gtools")#permutacie
library("mclust")

#Dunn index
Dunn<-function(x,cl, input){ #x su data, cl je clustrovaci vektor
  if(input!="vec" && input!="dif") stop("zadaj vec pre input s vektormi alebo dif pre maticu dolisnosti")
  
  X<-as.matrix(x)
  d<-dim(X)[2] #dim dat
  pocet<-dim(X)[1] #pocet dat
  k<-length(unique(cl))
  clust1<-matrix(0,nrow = pocet, ncol = pocet)
  clust2<-matrix(0,nrow = pocet, ncol = pocet)
  
  if(input=="vec") X<-as.matrix(dist(X, method = "euclidean")) #matica odlisnosti
  
  for (i in 1:pocet) {
    clust1[i,]<-as.numeric(cl!=cl[i]) #kontraindikatorova matica prislusnosti rovnakeho zhluku ako i-ty objekt
    clust2[i,]<-as.numeric(cl==cl[i]) #indikatorova matica prislusnosti rovnakeho zhluku ako i-ty objekt
  }
  
  X_a<-X*clust1
  X_b<-X*clust2
  X_a[X_a==0]<-NA
  X_b[X_b==0]<-NA
  a<-min(X_a,na.rm = TRUE)
  b<-max(X_b,na.rm = TRUE)
  
  
  s<-a/b
  return(s)
} 

#C-H index
C_H<-function(x,cl, centr=NULL){
  X<-as.matrix(x)
  d<-dim(X)[2] #dim dat
  pocet<-dim(X)[1] #pocet dat
  k<-length(unique(cl))
  centr<-matrix(0,nrow = k, ncol = d) #stred zhluku
  k_size<-rep(0,k)
  
  glob<-colMeans(X) #globalny stred
  b<-0
  for (i in 1:k) {
    centr[i,]<-colMeans(as.matrix(X[which(cl==i),]))
    k_size[i]<-length(cl[cl==i])
    
    stuff<-sqrt(diag(tcrossprod(t(t(as.matrix(X[which(cl==i),])) - centr[i,])))) #euklid. vzdialenost objektov od svojho stredu
    b<-b+sum(stuff)
  }
  a<-sum(sqrt(diag(tcrossprod((t(t(centr)-glob)))))*k_size)
  
  s<-((pocet-k)*a)/((k-1)*b)
  return(s)
}

#D-B index
D_B<-function(x,cl){
  X<-as.matrix(x)
  d<-dim(X)[2] #dim dat
  pocet<-dim(X)[1] #pocet dat
  k<-length(unique(cl))
  centr<-matrix(0,nrow = k, ncol = d)
  k_size<-rep(0,k)
  
  Err<-rep(0,k)
  for (i in 1:k) {
    centr[i,]<-colMeans(as.matrix(X[which(cl==i),]))
    k_size[i]<-length(cl[cl==i])
    
    stuff<-sqrt(diag(tcrossprod(t(t(as.matrix(X[which(cl==i),])) - centr[i,])))) #euklid. vzdialenost objektov od svojho stredu
    Err[i]<-sum(stuff)/k_size[i]
  }
  
  dif_centr<-as.matrix(dist(centr, method = "euclidean"))
  dif_centr[dif_centr==0]<-NA #pre istotu
  c<-rep(0,k)
  for (i in 1:k) {
    c[i]<-max((Err[i]+Err[-i])/dif_centr[i,-i])
  }
  
  s<-sum(c)/k
  return(s)
}

#statistika zhlukovanie, 1000 realizaci#############
#permutacie pre testovanie spravnosti zhlukovania
perm<-list(permutations(1,1),permutations(2,2),permutations(3,3), permutations(4,4), permutations(5,5), permutations(6,6), permutations(7,7), permutations(8,8), permutations(9,9))
perm_length<-(lengths(perm)/c(1,2,3,4,5,6,7,8,9))

n_gen<-sample(c(50:5000), 1000, replace = TRUE) #pocet datovych bodov v modelovych datach
k_gen<-sample(c(2:6), 1000, replace = TRUE)  #pocet zhlukov
d_gen<-sample(c(1:10), 1000, replace = TRUE) #dimenzia dat
while (sum(n_gen%%k_gen==0)!=length(n_gen)) { #delitelnost n a k pre jednoduchost, vychilenie pravdepodobnosti povazujeme za insignifikantne
  n_gen[n_gen%%k_gen!=0]<-n_gen[n_gen%%k_gen!=0]+1
}

#bench
bench_stat<-function(){
  #prazdne premenne
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(1, 10, by=0.1)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  #cyklus cez vsetky realizacie modelovych dat
  for (i in 1:1000) {
    #zhlukovanie
    dummy_kmeans<-kmeans(x_bench[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_bench[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(x_bench[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    
    for (h in 1:length(eps)) { #clustrovanie pre kazdy epsilon pri DBScan
      dummy_dbscan[h,]<-dbscan(x_bench[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    dummy_dbscan<-unique.matrix(dummy_dbscan) #unique zhlukovania DBScan nech nekontrolujeme to iste viackrat
    
    #priprava "spravneho" clustrovacieho vektoru (vektoru rozdelenia), musime ale pozriet vsetky permutacie mien zhlukov
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i]) #kazdy riadok ina permutacia mien zhlukov, tolko v riadku kolko je n
    for (j in 1:perm_length[k_gen[i]]) { 
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) ))) #sort lebo data su generovane postupne po zhlukoch a potrebujeme len permutacie zhlukov
    }
    
    #vyhodnotenie skore
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i] #najlepsie zhlukovanie z pohladu permutacie mien, percento spravne zaradenych datovych bodov
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) { 
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench)) #vektory uspesnosti zaradenia pre kazdu realizaciu dat
}

#noise
lownoise_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000) # _bench ostalo ako artefakt kopirovania, iny vyznam nema
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 6, by=0.1) #povodne zbehnute pre 1 az 10; neskor rozhodnute ze lepsie zacat od 0.1 ale z vysledkov sme uz vedeli ze staci len po 6; podobne aj pre dalsie funkcie
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) {
    dummy_x<-x_lownoise[[i]][[1]]
    n_new<-x_lownoise[[i]][[2]] #uprava kvoli zaokruhlovaniu z noise_gen()
    n_noise<-n_new-round(n_new*(1-(5/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise) #vektor velkosti zhlukov a noise
    
    dummy_kmeans<-kmeans(dummy_x, k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_new)
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(dummy_x, eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan) #urychlenie v pripade ze nepotrebujeme optim. epsilon; bolo vypocitane no nakoniec v praci nevyuzite
    
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new) #k+1 lebo aj noise, n_new lebo zaciatok funkcie
    for (j in 1:perm_length[k_gen[i]+1]) { 
      elong<-rep(perm[[k_gen[i]+1]][j,1]-1, n_weight[1]) #generovanie compare cez cyklus lebo nie kazdy zhluk ma rovnaky pocet datovych bodov
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g]-1, n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]-1) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_new
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)] #v praci nevyuzite
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

lownoise_stat_1<-function(){ #k+1 zadane, z tabuliek z prace ide o riadky s *
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  
  for (i in 1:1000) {  
    dummy_x<-x_lownoise[[i]][[1]]
    n_new<-x_lownoise[[i]][[2]]
    n_noise<-n_new-round(n_new*(1-(5/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise)
    dummy_kmeans<-kmeans(dummy_x, k_gen[i]+1, nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i]+1, diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=(k_gen[i]+1))$classification
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new)
    for (j in 1:perm_length[k_gen[i]+1]) {
      elong<-rep(perm[[k_gen[i]+1]][j,1], n_weight[1])
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g], n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench))
}

midnoise_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 6, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) { 
    dummy_x<-x_midnoise[[i]][[1]]
    n_new<-x_midnoise[[i]][[2]]
    n_noise<-n_new-round(n_new*(1-(15/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise)
    dummy_kmeans<-kmeans(dummy_x, k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_new)
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(dummy_x, eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new)
    for (j in 1:perm_length[k_gen[i]+1]) {
      elong<-rep(perm[[k_gen[i]+1]][j,1]-1, n_weight[1])
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g]-1, n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]-1) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_new
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps)) #list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench) 
}

midnoise_stat_1<-function(){ #k+1 zadane, z tabuliek z prace ide o riadky s *
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  
  for (i in 1:1000) { 
    dummy_x<-x_midnoise[[i]][[1]]
    n_new<-x_midnoise[[i]][[2]]
    n_noise<-n_new-round(n_new*(1-(15/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise)
    dummy_kmeans<-kmeans(dummy_x, k_gen[i]+1, nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i]+1, diss=FALSE, metric = "euclidean", cluster.only = TRUE) 
    dummy_mclust<-Mclust(dummy_x, G=(k_gen[i]+1))$classification
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new)
    for (j in 1:perm_length[k_gen[i]+1]) {
      elong<-rep(perm[[k_gen[i]+1]][j,1], n_weight[1])
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g], n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench))  
}

highnoise_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 10, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) {  
    dummy_x<-x_highnoise[[i]][[1]]
    n_new<-x_highnoise[[i]][[2]]
    n_noise<-n_new-round(n_new*(1-(35/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise)
    dummy_kmeans<-kmeans(dummy_x, k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_new)
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(dummy_x, eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new)
    for (j in 1:perm_length[k_gen[i]+1]) {
      elong<-rep(perm[[k_gen[i]+1]][j,1]-1, n_weight[1])
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g]-1, n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]-1) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_new
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

highnoise_stat_1<-function(){ #k+1 zadane, z tabuliek z prace ide o riadky s *
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  
  for (i in 1:1000) {  
    dummy_x<-x_highnoise[[i]][[1]]
    n_new<-x_highnoise[[i]][[2]]
    n_noise<-n_new-round(n_new*(1-(35/100)))
    n_weight<-c(rep((n_new-n_noise)/k_gen[i],k_gen[i]),n_noise)
    dummy_kmeans<-kmeans(dummy_x, k_gen[i]+1, nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i]+1, diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=(k_gen[i]+1))$classification
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]+1], ncol = n_new)
    for (j in 1:perm_length[k_gen[i]+1]) {
      elong<-rep(perm[[k_gen[i]+1]][j,1], n_weight[1])
      for (g in 2:(k_gen[i]+1)) {
        elong<-c(elong, rep(perm[[k_gen[i]+1]][j,g], n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]+1]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_new
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_new
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_new
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench)) 
}

#distr
distr_t_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 10, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) { 
    dummy_kmeans<-kmeans(x_distr_t[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_distr_t[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(x_distr_t[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(x_distr_t[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

distr_unif_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 3, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) { 
    dummy_kmeans<-kmeans(x_distr_unif[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_distr_unif[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(x_distr_unif[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(x_distr_unif[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

#covar
covar_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 10, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) {
    dummy_kmeans<-kmeans(x_covar[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_covar[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(x_covar[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(x_covar[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps)) #list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench) 
}

#near
near_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 7, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) {  
    dummy_kmeans<-kmeans(x_near[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_near[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE) 
    dummy_mclust<-Mclust(x_near[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(x_near[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

nearnear_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 5, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) { 
    dummy_kmeans<-kmeans(x_nearnear[[i]], k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(x_nearnear[[i]], k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(x_nearnear[[i]], G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(x_nearnear[[i]], eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(perm[[k_gen[i]]], nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      compares[j,]<-as.numeric(as.vector(sort( factor(compares[j,], levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps))
}

#weight
weight_stat<-function(){
  score_kmeans_bench<-rep.int(0, times = 1000)
  score_kmedoid_bench<-rep.int(0, times = 1000)
  score_mclust_bench<-rep.int(0, times = 1000)
  score_dbscan_bench<-rep.int(0, times = 1000)
  eps<-seq(0.1, 10, by=0.1)
  best_eps<-rep.int(0, times = 1000)
  maybe_score_dbscan_bench<-rep.int(0, times = length(eps))
  
  
  for (i in 1:1000) {
    dummy_x<-x_weight[[i]][[1]]
    n_weight<-x_weight[[i]][[2]]
    dummy_kmeans<-kmeans(dummy_x, k_gen[i], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster    
    dummy_kmedoid<-pam(dummy_x, k_gen[i], diss=FALSE, metric = "euclidean", cluster.only = TRUE)
    dummy_mclust<-Mclust(dummy_x, G=k_gen[i])$classification
    dummy_dbscan<-matrix(0, nrow = length(eps), ncol = n_gen[i])
    for (h in 1:length(eps)) {
      dummy_dbscan[h,]<-dbscan(dummy_x, eps = eps[h], minPts = 2*d_gen[i])$cluster
    }
    #dummy_dbscan<-unique.matrix(dummy_dbscan)
    
    compares<-matrix(0, nrow = perm_length[k_gen[i]], ncol = n_gen[i])
    for (j in 1:perm_length[k_gen[i]]) {
      elong<-rep(perm[[k_gen[i]]][j,1], n_weight[1]) #vid. lownoise pre vysvetlenie
      for (g in 2:(k_gen[i])) {
        elong<-c(elong, rep(perm[[k_gen[i]]][j,g], n_weight[g]))
      }
      compares[j,]<-as.numeric(as.vector(sort( factor(elong, levels = perm[[k_gen[i]]][j,]) )))
    }
    
    score_kmeans_bench[i]<-max(colSums(t(compares)==dummy_kmeans))/n_gen[i]
    score_kmedoid_bench[i]<-max(colSums(t(compares)==dummy_kmedoid))/n_gen[i]
    score_mclust_bench[i]<-max(colSums(t(compares)==dummy_mclust))/n_gen[i]
    for (g in 1:dim(dummy_dbscan)[1]) {
      maybe_score_dbscan_bench[g]<-max(colSums(t(compares)==dummy_dbscan[g,]))/n_gen[i]
    }
    score_dbscan_bench[i]<-max(maybe_score_dbscan_bench)
    best_eps[i]<-eps[which.max(maybe_score_dbscan_bench)]
  }
  return(list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench, best_eps)) #list(score_kmeans_bench, score_kmedoid_bench, score_mclust_bench, score_dbscan_bench) 
}



#statistika urcovania poctu dat, 1000 realizaci#########

#prepisanie zhlukovacich funkcii pre Gap statistic, vyzaduje prvy input data, duhy input pocet zhlukov, output vektor clustrovania
kmeans1 <- function(x,k) list(cluster = kmeans(x,k, nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster)
pam1 <- function(x,k) list(cluster = pam(x,k, diss=FALSE, metric = "euclidean", cluster.only=TRUE))
Mclust1 <- function(x,k) list(cluster = Mclust(x,G=k)$classification)

#bench
bench_urcenie<-function(){
  k_now<-c(2:6) #cez ake pocty zhlukov idem hladat
  #priprava prazdnych premennych
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  #cyklus cez vsetky realizacie modelovych dat
  for (i in 1:1000) {
    #Gap statistic pozrie automaticky cez vsetky k od 1 po K.max
    dummy_gap<-clusGap(x_bench[[i]], FUNcluster = kmeans1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    #zhlukovanie cez vsetky k
    for (j in 2:6) {
      #MClust, ulozim aic, bic
      dummy_mclust<-Mclust(x_bench[[i]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){ #velmi ojedinele error "replacment has length zero", jedna z podmienok v if to vyriesila. Neviem preco ani ktora
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic #Mclust pouziva -BIC takze treba napravit znamienko
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik
      }
      
      #vytazne zhlukovanie, nasledne ulozenie hodnoty indexov
      dummy_best_clust<-kmeans(x_bench[[i]], k_now[k_now==j], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster
      dunn[j-1]<-Dunn(x_bench[[i]], dummy_best_clust, "vec")
      ch[j-1]<-C_H(x_bench[[i]], dummy_best_clust) 
      db[j-1]<-D_B(x_bench[[i]], dummy_best_clust) 
      sil[j-1]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_bench[[i]]))))$avg.width
    }
    
    #odhad na zaklade hodnot ukazovatelov; +1 vsade lebo k ide od 2 a teda na prvom mieste vektora bic je hodnota pre 2 zhluky. which.min teda vrati k-1 pre optim k
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap)) #vektory odhadnuteho poctu zhlukov
}

#priklad vypocitania uspesnosti BIC pri benchmarku v percentach
test<-bench_urcenie()
mean(test[[1]]==k_gen)*100


#noise
lownoise_urcenie<-function(){
  k_now<-c(2:6) 
  eps_now<-seq(from = 0.1, to = 6, by = 0.1)
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(eps_now))
  ch<-rep(NaN, length(eps_now))
  db<-rep(NaN, length(eps_now))
  sil<-rep(NaN, length(eps_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_lownoise[[i]][[1]], FUNcluster = pam1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_lownoise[[i]][[1]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik  
      }
    }
    
    #verzia bez sumu, v tabulkach z prace ide teda o verziu BEZ *
    for (j in 1:length(eps_now)) { 
      dummy_best_clust<-dbscan(x_lownoise[[i]][[1]], eps_now[j], 2*d_gen[i])$cluster
      X<-x_lownoise[[i]][[1]][dummy_best_clust!=0,] #iba data ktore niesu noise
      dummy_best_clust<-dummy_best_clust[dummy_best_clust!=0] #clustrovaci vektor bez noise
      if(length(unique(dummy_best_clust))>1){ #indexi definovane len pre 2+ zhluky
        dunn[j]<-Dunn(X, dummy_best_clust, "vec")
        ch[j]<-C_H(X, dummy_best_clust) 
        db[j]<-D_B(X, dummy_best_clust) 
        sil[j]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(X))))$avg.width
      }
    }
    
    #verzia so sumom, v tabulkach z prace ide teda o verziu S *
    #for (j in 1:length(eps_now)) { 
    #  dummy_best_clust<-dbscan(x_lownoise[[i]][[1]], eps_now[j], 2*d_gen[i])$cluster
    #  if(0%in%dummy_best_clust) dummy_best_clust<-dummy_best_clust+1 #ak je noise tak pripocitame k menam zhlukov +1, indexi nevedia pocitat so zhlukom s nazvom<1
    #  if(length(unique(dummy_best_clust))>1){          #indexi definovane len pre 2+ zhluky
    #    dunn[j]<-Dunn(x_lownoise[[i]][[1]], dummy_best_clust, "vec")
    #    ch[j]<-C_H(x_lownoise[[i]][[1]], dummy_best_clust) 
    #    db[j]<-D_B(x_lownoise[[i]][[1]], dummy_best_clust) 
    #    sil[j]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_lownoise[[i]][[1]]))))$avg.width
    #  }
    #}
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-eps_now[which.max(dunn)] #iba odhad optimalneho epsilon, pcoet zhlukov treba najst potom; vid nizsie
    guess_ch[i]<-eps_now[which.max(ch)]
    guess_db[i]<-eps_now[which.min(db)]
    guess_sil[i]<-eps_now[which.max(sil)]
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}

#priklad zistenia optimalneho poctu zhlukov z optimalneho parametra epsilon
stat_urcenie_lownoise_index_nonoise<-list(rep(NA, 1000), rep(NA, 1000), rep(NA, 1000), rep(NA, 1000))
for (i in 1:1000) {
  k<-dbscan(x_lownoise[[i]][[1]], hope_urcenie_lownoise_index_nonoise_eps[[1]][i])$cluster
  stat_urcenie_lownoise_index_nonoise[[1]][i]<-length(unique(k[k!=0]))
  
  k<-dbscan(x_lownoise[[i]][[1]], hope_urcenie_lownoise_index_nonoise_eps[[2]][i])$cluster
  stat_urcenie_lownoise_index_nonoise[[2]][i]<-length(unique(k[k!=0]))
  
  k<-dbscan(x_lownoise[[i]][[1]], hope_urcenie_lownoise_index_nonoise_eps[[3]][i])$cluster
  stat_urcenie_lownoise_index_nonoise[[3]][i]<-length(unique(k[k!=0]))
  
  k<-dbscan(x_lownoise[[i]][[1]], hope_urcenie_lownoise_index_nonoise_eps[[4]][i])$cluster
  stat_urcenie_lownoise_index_nonoise[[4]][i]<-length(unique(k[k!=0]))
}


midnoise_urcenie<-function(){
  k_now<-c(2:6) 
  eps_now<-seq(from = 0.1, to = 6, by = 0.1) 
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(eps_now))
  ch<-rep(NaN, length(eps_now))
  db<-rep(NaN, length(eps_now))
  sil<-rep(NaN, length(eps_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_midnoise[[i]][[1]], FUNcluster = Mclust1, K.max = 7, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_midnoise[[i]][[1]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik  
      }
    }
    
    
    #verzia bez sumu, v tabulkach z prace ide teda o verziu BEZ *
    for (j in 1:length(eps_now)) { 
      dummy_best_clust<-dbscan(x_midnoise[[i]][[1]], eps_now[j], 2*d_gen[i])$cluster
      X<-x_midnoise[[i]][[1]][dummy_best_clust!=0,] #iba data ktore niesu noise
      dummy_best_clust<-dummy_best_clust[dummy_best_clust!=0] #clustrovaci vektor bez noise
      if(length(unique(dummy_best_clust))>1){ #indexi definovane len pre 2+ zhluky
        dunn[j]<-Dunn(X, dummy_best_clust, "vec")
        ch[j]<-C_H(X, dummy_best_clust) 
        db[j]<-D_B(X, dummy_best_clust) 
        sil[j]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(X))))$avg.width
      }
    }
    
    #pre verziu so sumom vid. lownoise
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-eps_now[which.max(dunn)]
    guess_ch[i]<-eps_now[which.max(ch)]
    guess_db[i]<-eps_now[which.min(db)]
    guess_sil[i]<-eps_now[which.max(sil)]
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}


highnoise_urcenie<-function(){
  k_now<-c(2:6) 
  eps_now<-seq(from = 0.1, to = 10, by = 0.1)
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(eps_now))
  ch<-rep(NaN, length(eps_now))
  db<-rep(NaN, length(eps_now))
  sil<-rep(NaN, length(eps_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_highnoise[[i]][[1]], FUNcluster = pam1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_highnoise[[i]][[1]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik 
      }
    }
    
    
    #verzia bez sumu, v tabulkach z prace ide teda o verziu BEZ *
    for (j in 1:length(eps_now)) { 
      dummy_best_clust<-dbscan(x_highnoise[[i]][[1]], eps_now[j], 2*d_gen[i])$cluster
      X<-x_highnoise[[i]][[1]][dummy_best_clust!=0,] #iba data ktore niesu noise
      dummy_best_clust<-dummy_best_clust[dummy_best_clust!=0] #clustrovaci vektor bez noise
      if(length(unique(dummy_best_clust))>1){ #indexi definovane len pre 2+ zhluky
        dunn[j]<-Dunn(X, dummy_best_clust, "vec")
        ch[j]<-C_H(X, dummy_best_clust) 
        db[j]<-D_B(X, dummy_best_clust) 
        sil[j]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(X))))$avg.width
      }
    }
    
    #pre verziu so sumom vid. lownoise
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-eps_now[which.max(dunn)]
    guess_ch[i]<-eps_now[which.max(ch)]
    guess_db[i]<-eps_now[which.min(db)]
    guess_sil[i]<-eps_now[which.max(sil)]
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}


#distr
distr_t_urcenie<-function(){
  k_now<-c(2:6)
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_distr_t[[i]], FUNcluster = pam1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_distr_t[[i]], G=k_now[k_now==j], verbose = FALSE)
      #if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){ #tu z nejakeho dovodu nebolo treba
      #  bic[j-1]<-NaN
      #  aic[j-1]<-NaN
      #}
      #else{
      bic[j-1]<-(-1)*dummy_mclust$bic
      aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik 
      #}
      
      dummy_best_clust<-pam(x_distr_t[[i]], k_now[k_now==j],diss=FALSE, metric = "euclidean", cluster.only = TRUE)
      dunn[j-1]<-Dunn(x_distr_t[[i]], dummy_best_clust, "vec")
      ch[j-1]<-C_H(x_distr_t[[i]], dummy_best_clust) 
      db[j-1]<-D_B(x_distr_t[[i]], dummy_best_clust) 
      sil[j-1]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_distr_t[[i]]))))$avg.width
    }
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}

distr_unif_urcenie<-function(){
  k_now<-c(2:6) 
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_distr_unif[[i]], FUNcluster = kmeans1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_distr_unif[[i]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik
      }
      
      dummy_best_clust<-kmeans(x_distr_unif[[i]], k_now[k_now==j], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster
      dunn[j-1]<-Dunn(x_distr_unif[[i]], dummy_best_clust, "vec")
      ch[j-1]<-C_H(x_distr_unif[[i]], dummy_best_clust) 
      db[j-1]<-D_B(x_distr_unif[[i]], dummy_best_clust) 
      sil[j-1]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_distr_unif[[i]]))))$avg.width
    }
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}

#covar
covar_urcenie<-function(){
  k_now<-c(2:6) 
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_covar[[i]], FUNcluster = Mclust1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_covar[[i]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik
      }
      
      dunn[j-1]<-Dunn(x_covar[[i]], dummy_mclust$classification, "vec")
      ch[j-1]<-C_H(x_covar[[i]], dummy_mclust$classification) 
      db[j-1]<-D_B(x_covar[[i]], dummy_mclust$classification) 
      sil[j-1]<-summary(silhouette(dummy_mclust$classification, dmatrix = as.matrix(dist(x_covar[[i]]))))$avg.width
    }
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}


#near/nearnear
near_urcenie<-function(){
  k_now<-c(2:6) 
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_near[[i]], FUNcluster = kmeans1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_near[[i]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik  
      }
      
      dummy_best_clust<-kmeans(x_near[[i]], k_now[k_now==j], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster
      dunn[j-1]<-Dunn(x_near[[i]], dummy_best_clust, "vec") 
      ch[j-1]<-C_H(x_near[[i]], dummy_best_clust) 
      db[j-1]<-D_B(x_near[[i]], dummy_best_clust) 
      sil[j-1]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_near[[i]]))))$avg.width
    }
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}

nearnear_urcenie<-function(){
  k_now<-c(2:6) 
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(k_now))
  ch<-rep(NaN, length(k_now))
  db<-rep(NaN, length(k_now))
  sil<-rep(NaN, length(k_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_nearnear[[i]], FUNcluster = kmeans1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_nearnear[[i]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik 
      }
      
      dummy_best_clust<-kmeans(x_nearnear[[i]], k_now[k_now==j], nstart = 50, iter.max = 50, algorithm = "Lloyd")$cluster
      dunn[j-1]<-Dunn(x_nearnear[[i]], dummy_best_clust, "vec") 
      ch[j-1]<-C_H(x_nearnear[[i]], dummy_best_clust) 
      db[j-1]<-D_B(x_nearnear[[i]], dummy_best_clust) 
      sil[j-1]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_nearnear[[i]]))))$avg.width
    }
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-which.max(dunn)+1
    guess_ch[i]<-which.max(ch)+1
    guess_db[i]<-which.min(db)+1
    guess_sil[i]<-which.max(sil)+1
  }
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}

#weight
weight_urcenie<-function(){
  k_now<-c(2:6) 
  eps_now<-seq(from = 0.1, to = 10, by = 0.1)
  bic<-rep(NaN, length(k_now))
  aic<-rep(NaN, length(k_now))
  dunn<-rep(NaN, length(eps_now))
  ch<-rep(NaN, length(eps_now))
  db<-rep(NaN, length(eps_now))
  sil<-rep(NaN, length(eps_now))
  guess_gap<-rep(NaN,1000)
  guess_bic<-rep(NaN,1000)
  guess_aic<-rep(NaN,1000)
  guess_dunn<-rep(NaN,1000)
  guess_ch<-rep(NaN,1000)
  guess_db<-rep(NaN,1000)
  guess_sil<-rep(NaN,1000)
  
  for (i in 1:1000) {
    dummy_gap<-clusGap(x_weight[[i]][[1]], FUNcluster = kmeans1, K.max = 6, d.power = 2, verbose = FALSE, B=50)$Tab
    guess_gap[i]<-maxSE(dummy_gap[,3],dummy_gap[,4], method = "Tibs2001SEmax")
    
    
    for (j in 2:6) {
      dummy_mclust<-Mclust(x_weight[[i]][[1]], G=k_now[k_now==j], verbose = FALSE)
      if(length(dummy_mclust$bic)==0 || is.numeric(dummy_mclust$bic)==FALSE || is.na(dummy_mclust$bic)){
        bic[j-1]<-NaN
        aic[j-1]<-NaN
      }
      else{
        bic[j-1]<-(-1)*dummy_mclust$bic
        aic[j-1]<-2*dummy_mclust$df-2*dummy_mclust$loglik 
      }
    }
    
    #verzia pre noise
    for (j in 1:length(eps_now)) { 
      dummy_best_clust<-(dbscan(x_weight[[i]][[1]], eps_now[j], 2*d_gen[i])$cluster)+1
      if(length(unique(dummy_best_clust))!=1){
        dunn[j]<-Dunn(x_weight[[i]][[1]], dummy_best_clust, "vec")
        ch[j]<-C_H(x_weight[[i]][[1]], dummy_best_clust) 
        db[j]<-D_B(x_weight[[i]][[1]], dummy_best_clust) 
        sil[j]<-summary(silhouette(dummy_best_clust, dmatrix = as.matrix(dist(x_weight[[i]][[1]]))))$avg.width
      }
    }
    
    #pre verziu bez noise vid. lownoise funkciu
    
    guess_bic[i]<-which.min(bic)+1
    guess_aic[i]<-which.min(aic)+1
    guess_dunn[i]<-eps_now[which.max(dunn)]
    guess_ch[i]<-eps_now[which.max(ch)]
    guess_db[i]<-eps_now[which.min(db)]
    guess_sil[i]<-eps_now[which.max(sil)]
  }
  
  return(list(guess_bic,guess_aic,guess_dunn,guess_ch,guess_db,guess_sil,guess_gap))
}
