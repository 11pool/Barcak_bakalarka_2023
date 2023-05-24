#Kód v jazyku R pouitı pre bakalársku prácu:
#Urèovanie poètu zhlukov v dátach
#Adam Barèák
#Univerzita Komenského FMFI, 2023
#####################################

bench_gen<-function(n,k,d){ #n musi byt delitelne k, vektory v stlpcoch ale return dava v riadkoch
  check<-0 #check ci su stredy zhlukov dobre vygenerovane
  x<-matrix(0,nrow = d, ncol = n)
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4) # heurisitcky zvolene polstrany kocky z ktorej generujeme stredy pre rozne dimenzie
  while (check==0) {
    stred<-matrix(runif(k*d, -range[d], range[d]), nrow = d, ncol = k)  #generovanie stredov zhlukov
    if(all(dist(t(stred))>10)){   #min vzdialenos stredov; heuristicky zvolena hodnota 10 ako vzdialenost ked sa este medzi zhluky zmesti iny
      check<-1
    }
  }
  
  if(k!=1){ #nakoniec k=1 nepouzite
    for (j in 0:(k-1)) {
      x[,(j*n/k+1):((j+1)*n/k)]<-matrix(rnorm((n/k)*d, 0, 1),nrow=d)+stred[,j+1]  #covar. matica zvolena ako identitia; vygenerujem a posuniem do spravneho stredu
    }
  }
  return(t(x))
}

noise_gen<-function(n,k,d, frac){ #vektory v stlpcoch ale return dava v riadkoch
  check<-0
  while (round(n*(1-frac))%%k!=0 || n%%(k+1)!=0) { #minimalne navysenie poctu datovych bodov aby nebol problem s delenim. Prva podmienka je aby sa dali zvysne data mimo sum rozdelit rovnomerne. Druha osetruje (k+1) statistiku
    n<-n+1
  }
  n1<-round(n*(1-frac))
  x<-matrix(0,nrow = d, ncol = n)
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4)
  while (check==0) {
    stred<-matrix(runif(k*d, -range[d], range[d]), nrow = d, ncol = k)  
    if(all(dist(t(stred))>10)){ 
      check<-1
    }
  }
  
  for (j in 0:(k-1)) {
    x[,(j*n1/k+1):((j+1)*n1/k)]<-matrix(rnorm((n1/k)*d, 0, 1),nrow=d)+stred[,j+1]
  }
  
  x[,((n1)+1):n]<-runif((n-(n1))*d, -range[d]-10, range[d]+10) #+10 nech zahrnie pripad ze zhluky vytrcaju z kocky s polhranou range
  return(list(t(x),n))
}

distr_gen<-function(n,k,d, distr){ #n delitelne k, vektory v stlpcoch ale return dava v riadkoch; distr bud "t" alebo "unif" alebo "norm"(co je ale bench)
  check<-0
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4)
  min_dist<-10 #rovnaka hodnota len parametrizovane
  if(distr=="unif") { #bolo treba zmensit lebo velmi daleko od seba boli pre uniform
    min_dist<-3
    range<-c(10, 4, 4, 2, 2, 2, 2, 2, 2, 2)
  }
  while (check==0) {
    stred<-matrix(runif(k*d, -range[d], range[d]), nrow = d, ncol = k)
    if(all(dist(t(stred))>min_dist)){ 
      check<-1
    }
  }
  
  if(distr=="t"){   #t-distr
    x<-matrix(0,nrow = d, ncol = n)
    u<-rchisq(n,3) #chi^2, deg=3
    u<-sqrt(3/u) #deg=3
    for (j in 0:(k-1)) {
      x[,(j*n/k+1):((j+1)*n/k)]<-t(t(matrix(rnorm((n/k)*d, 0, 1),nrow=d))*u[(j*n/k+1):((j+1)*n/k)])+stred[,j+1] #(norm a chi^2)+posun -> t-distr so stredom v posun
    }
  }
  
  if(distr=="norm"){   #norm-distr
    x<-matrix(0,nrow = d, ncol = n)
    for (j in 0:(k-1)) {
      x[,(j*n/k+1):((j+1)*n/k)]<-matrix(rnorm((n/k)*d, 0, 1),nrow=d)+stred[,j+1]
    }
  }
  
  if(distr=="unif"){   #unif-distr; princip generovania popisany v praci
    x<-matrix(rnorm(n*(d+2), 0, 1),nrow = d+2, ncol = n)
    size<-sqrt(diag(crossprod(x)))
    x<-t(t(x)/size)
    x<-x[1:d,] #toto je rovnomerne na guli
    
    if(is.vector(x)==TRUE) x<-t(as.matrix(x)) #v pripade d=1 aby vsetko sedelo
    
    for (j in 0:(k-1)) {
      x[,(j*n/k+1):((j+1)*n/k)]<-x[,(j*n/k+1):((j+1)*n/k)]+stred[,j+1] #rovnomerne v guli + posun
    }
  }
  
  return(t(x))
}

covar_gen<-function(n,k,d){
  check<-0
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4)
  x<-matrix(0,nrow = d, ncol = n)
  while (check==0) {
    stred<-matrix(runif(k*d, -range[d], range[d]), nrow = d, ncol = k)
    if(all(dist(t(stred))>10)){ 
      check<-1
    }
  }
  
  if(d==1){
    a<-runif(k, 0.1, 10) #potencialne disperzie
    for (j in 0:(k-1)) {
      x[,(j*n/k+1):((j+1)*n/k)]<-matrix(rnorm((n/k)*d, 0, a[j+1]),nrow=d)+stred[,j+1]
    }
  }
  
  else{
    for (j in 0:(k-1)) {
      eigen_values<-diag(sort(runif(d, 0.1, 10), decreasing = TRUE)) #zoradene vlastne hodnoty z rovnakeho intervalu ako pre d=1
      eigen_vectors<-qr.Q(qr(matrix(rnorm(d^2), d))) #vlastn. vektory rovnomerne na guli(a teda rovnomerne do kazdej strany ukazujuce)
      A<-t(eigen_vectors)%*%eigen_values%*%eigen_vectors #random covar matica
      x[,(j*n/k+1):((j+1)*n/k)]<-t(t(chol(A))%*%matrix(rnorm((n/k)*d, 0, 1),nrow=d))+stred[,j+1]
    }
  }
  return(t(x))
}

nearnear_gen<-function(n,k,d, min_dist, max_dist, in_stred){
  check<-0
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4) 
  x<-matrix(0,nrow = d, ncol = n)
  stred<-matrix(0, nrow = d, ncol = k)
  while (check==0) {
    stred[,1:(k-in_stred)]<-matrix(runif((k-in_stred)*d, -range[d], range[d]), nrow = d, ncol = (k-in_stred))  #normalne generovanie k-in_stred zhlukov
    if(all(dist(t(stred)[1:(k-in_stred),])>min_dist) || (k-in_stred)==1){
      check<-1
    }
  }
  for (i in 1:in_stred) { #in_stred zhluky sa daju blizko uz vygenerovanym stredom
    a<-runif(d)
    stred[,(k-in_stred+i)]<-stred[,runif(1,1,k-in_stred)]+(a/length(a))*max_dist
  }
  
  
  for (j in 0:(k-1)) {
    x[,(j*n/k+1):((j+1)*n/k)]<-matrix(rnorm((n/k)*d, 0, 1),nrow=d)+stred[,j+1]
  }
  return(t(x))
}

near_gen<-function(n,k,d, min_dist){  
  check<-0
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4) 
  x<-matrix(0,nrow = d, ncol = n)
  stred<-matrix(0, nrow = d, ncol = k)
  while (check==0) {
    stred[,1:(k)]<-matrix(runif((k)*d, -range[d], range[d]), nrow = d, ncol = (k)) 
    if(all(dist(t(stred))>min_dist)){
      check<-1
    }
  }
  
  for (j in 0:(k-1)) {
    x[,(j*n/k+1):((j+1)*n/k)]<-matrix(rnorm((n/k)*d, 0, 1),nrow=d)+stred[,j+1]  
  }
  return(t(x))
}

weight_gen<-function(n,k,d){ 
  check<-0
  range<-c(35, 15, 10, 8, 7, 6, 5, 5, 5, 4) 
  x<-matrix(0,nrow = d, ncol = n)
  while (check==0) {
    stred<-matrix(runif(k*d, -range[d], range[d]), nrow = d, ncol = k)  
    if(all(dist(t(stred))>10)){
      check<-1
    }
  }
  
  n_weight<-rep(0,k) 
  while(sum(n_weight)!=n || sum(n_weight==0)!=0){ #pre istotu
    weight<-runif(k,0,1)  #rovnomerne vygenerovane vahy zhlukov
    weight<-weight/sum(weight)
    n_weight<-round(n*weight) #norm. na sucet vah = 1 a zaokruhlenie lebo celociselne pocty
  }
  
  
  x[,1:(n_weight[1])]<-matrix(rnorm((n_weight[1])*d, 0, 1),nrow=d)+stred[,1]
  for (j in 1:(k-1)) {
    x[,(sum(n_weight[1:(j)])+1):(sum(n_weight[1:(j+1)]))]<-matrix(rnorm((n_weight[j+1])*d, 0, 1),nrow=d)+stred[,j+1]
  }
  return(list(t(x), n_weight))
}

#priklad generovania dat pouzitim funkcie
x_weight<-vector(mode = "list", length = 1000)
for(i in 1:1000){
  x_weight[[i]]<-weight_gen(n_gen[i], k_gen[i], d_gen[i]) 
}