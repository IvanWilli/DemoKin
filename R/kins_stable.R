# base in Caswell (2019)
### main function
# it needs:
  # p: a survival probability vector. Its lenght is the number of age stages
  # f: a fertility vector. Same length than p.
  # age: simple ages by default
  # cum_deaths: tell if cumulate death experience

kins_stable <- function(p, f, age = 0:(length(p)-1),
                       cum_deaths = F, birth_female = 1/2.04){

  # make matrix transition from vectors
  ages = length(age)
  Ut = Mt = zeros = Dcum = matrix(0, nrow=ages, ncol=ages)
  Ut[row(Ut)-1 == col(Ut)] <- p[-ages]
  Ut[ages, ages] = p[ages]
  diag(Mt) = 1 - p
  if(cum_deaths){
    diag(Dcum) = 1
  }
  Ut = rbind(cbind(Ut,zeros),
             cbind(Mt,Dcum))
  Ft = matrix(0, nrow=ages*2, ncol=ages*2)

  # Caswell's assumption
  Ft[1,1:ages] = f * p * birth_female

  # stable age distr
  A = Ut[1:ages,1:ages] + Ft[1:ages,1:ages]
  A_decomp = eigen(A)
  lambda = as.double(A_decomp$values[1])
  w = as.double(A_decomp$vectors[,1])/sum(as.double(A_decomp$vectors[,1]))
  pi = w*A[1,]/sum(w*A[1,])

  # identity
  e = matrix(0, ages * 2, ages * 2)
  diag(e[1:ages,1:ages]) = 1

  # initial vectors
  d = gd = m = gm = ggm = os = ys = nos = nys = oa = ya = coa = cya = matrix(0, ages * 2, ages)

  m[,1] = c(pi, rep(0,ages))
  for(i in 1:(ages-1)){
    d[,i+1]   = Ut %*% d[,i] + Ft %*% e[,i]
    gd[,i+1]  = Ut %*% gd[,i] + Ft %*% d[,i]
    m[,i+1]   = Ut %*% m[,i]
    ys[,i+1]  = Ut %*% ys[,i] + Ft %*% m[,i]
    nys[,i+1] = Ut %*% nys[,i] + Ft %*% ys[,i]
  }

  gm[1:ages,1] = m[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    gm[,i+1]  = Ut %*% gm[,i]
  }

  ggm[,1] = gm[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    ggm[,i+1]  = Ut %*% ggm[,i]
  }

  os[,1]  = d[1:ages,] %*% pi
  nos[,1] = gd[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    os[,i+1]  = Ut %*% os[,i]
    nos[,i+1] = Ut %*% nos[,i] + Ft %*% os[,i]
  }

  oa[,1]  = os[1:ages,] %*% pi
  ya[,1]  = ys[1:ages,] %*% pi
  coa[,1] = nos[1:ages,] %*% pi
  cya[,1] = nys[1:ages,] %*% pi
  for(i in 1:(ages-1)){
    oa[,i+1]  = Ut %*% oa[,i]
    ya[,i+1]  = Ut %*% ya[,i]  + Ft %*% gm[,i]
    coa[,i+1] = Ut %*% coa[,i] + Ft %*% oa[,i]
    cya[,i+1] = Ut %*% cya[,i] + Ft %*% ya[,i]
  }

  # get results
  kins_list <- list(d=d,gd=gd,m=m,gm=gm,ggm=ggm,os=os,ys=ys,
                    nos=nos,nys=nys,oa=oa,ya=ya,coa=coa,cya=cya)
  kins <- map2(kins_list, names(kins_list),
               function(x,y){
                    out = as.data.frame(x)
                    colnames(out) = age
                    out %>%
                      mutate(kin = y,
                             age_kin = rep(age,2),
                             alive = c(rep("yes",ages), rep("no",ages))) %>%
                      gather(age_ego,count,-age_kin, -kin, -alive) %>%
                      mutate(age_ego = as.integer(age_ego)) %>%
                      rename(x = age_ego,
                             x_kin = age_kin)
                    }
               ) %>%
              reduce(rbind)

  # return or not death experience
  if(cum_deaths == FALSE){
    kins <- kins %>%
              filter(alive == "yes") %>%
              select(-alive)
  }

  return(kins)
}
