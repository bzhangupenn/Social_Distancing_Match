# Programs for strengthening IV project, applying it to data
# Code for matching for Donut Hole IV project

###############################################################################################
# Functions to be used
library(RItools)
library(nbpMatching)
# rank based Mahalanobis distance between each pair
smahal=function(X){
  X<-as.matrix(X)
  n<-dim(X)[1]
  k<-dim(X)[2]
  for (j in 1:k) X[,j]<-rank(X[,j])
  cv<-cov(X)
  vuntied<-var(1:n)
  rat<-sqrt(vuntied/diag(cv))
  cv<-diag(rat)%*%cv%*%diag(rat)
  out<-matrix(NA,n,n)
  library(MASS)
  icov<-ginv(cv)
  for (i in 1:n) out[i,]<-mahalanobis(X,X[i,],icov,inverted=T)
  out
}

#Wrap function for matching function.  Creates sinks.
nbp_match_reverse_caliper <- function(select1row, X, IV, mindist, sinks){

  #Indices of columns with variation
  cols.w.var = which(apply(X[select1row,],2,var) > 0)

  #Create distance matrix
  distmat = smahal(X[select1row, cols.w.var]);
  distmat=round(distmat);

  #Make sure we don't match babies who have similar encouragement.
  distmat.adjust <- min(sd(distmat))
  dtt <- IV[select1row]
  store <- matrix(dtt, nrow = length(dtt), ncol=length(dtt))
  store <- abs(store - t(store))
  distmat[store < mindist] = distmat[store < mindist] + ifelse((max(store[store<mindist])-store[store<mindist])>2,
                             2,(max(store[store<mindist])-store[store<mindist]))*distmat.adjust*10


  #Encourage matching clusters with similar percen of male
  income = X[select1row, 'median_house_income']
  store_inc <- matrix(income, nrow = length(income), ncol=length(income))
  store_inc <- abs(store_inc - t(store_inc))
  distmat[store_inc > 800] = distmat[store_inc > 800] + 10

  black = X[select1row, 'bl_hisp']
  store_bl <- matrix(black, nrow = length(black), ncol=length(black))
  store_bl <- abs(store_bl - t(store_bl))
  distmat[store_bl > 0.1] = distmat[store_bl > 0.1] + 10

  female = X[select1row, 'female']
  store_fem <- matrix(female, nrow = length(female), ncol=length(female))
  store_fem <- abs(store_fem - t(store_fem))
  distmat[store_fem > 0.05] = distmat[store_fem > 0.05] + 8

  above_65 = X[select1row, 'above_65']
  store_65 <- matrix(above_65, nrow = length(above_65), ncol=length(above_65))
  store_65 <- abs(store_65 - t(store_65))
  distmat[store_65 > 0.05] = distmat[store_65 > 0.05] + 7

  #
  #Encourage matching clusters with similar percen of male
  traffic = X[select1row, 'traffic']
  store_traf <- matrix(traffic, nrow = length(traffic), ncol=length(traffic))
  store_traf <- abs(store_traf - t(store_traf))
  distmat[store_traf > 30] = distmat[store_traf > 30] + 20
  #
  #Encourage matching clusters with similar percen of male
  college = X[select1row, 'some_college']
  store_col <- matrix(college, nrow = length(college), ncol=length(college))
  store_col <- abs(store_col - t(store_col))
  distmat[store_col > 0.03] = distmat[store_col > 0.03] + 22

  smoke = X[select1row, 'smoking']
  store_smoke <- matrix(smoke, nrow = length(smoke), ncol=length(smoke))
  store_smoke <- abs(store_smoke - t(store_smoke))
  distmat[store_smoke > 0.05] = distmat[store_smoke > 0.05] + 5

  rucc = X[select1row, 'RUCC']
  store_rucc <- matrix(rucc, nrow = length(rucc), ncol=length(rucc))
  store_rucc <- abs(store_rucc - t(store_rucc))
  distmat[store_rucc > 2] = distmat[store_rucc > 2] + 10
  # #
  # #Encourage matching clusters with similar percen of male
  # exp_1 = X[select1row, 'avg_exp_ill_0316_0322']
  # store_exp1 <- matrix(exp_1, nrow = length(exp_1), ncol=length(exp_1))
  # store_exp1 <- abs(store_exp1 - t(store_exp1))
  # distmat[store_exp1 > 0.1] = distmat[store_exp1 > 0.1] + 10

  #Encourage matching clusters with similar percen of male
  exp_2 = X[select1row, 'avg_exp_ill_0330_0405']
  store_exp2 <- matrix(exp_2, nrow = length(exp_2), ncol=length(exp_2))
  store_exp2 <- abs(store_exp2 - t(store_exp2))
  distmat[store_exp2 > 0.007] = distmat[store_exp2 > 0.007] + 30
  # #
  #Encourage matching clusters with similar percen of male
  density = X[select1row, 'population.density']
  store_den <- matrix(density, nrow = length(density), ncol=length(density))
  store_den <- abs(store_den - t(store_den))
  distmat[store_den > 200] = distmat[store_den > 200] + 26




  #Create artifical sinks
  size <- dim(distmat)[2]
  num.sinks <- size*sinks
  num.sinks <- 2*ceiling((size+num.sinks)/2) - size
  total <- size + num.sinks
  distmat <- cbind(rbind(distmat,matrix(0,nrow=num.sinks,ncol=size)),matrix(0,nrow=total,ncol=num.sinks));
  if(num.sinks>0){distmat[(size+1):total,(size+1):total] <- max(distmat)*3}


  #Do the matching
  distmat=round(distmat);
  distmat=distancematrix(distmat)
  matching=nonbimatch(distmat)$matches;



  return(matching)
}

# This function takes in a matching object returned
# from nbpmatching, and return a list of three
# objects: encouraged dataset, control dataset

#####################################################################################################
###The data for matching

dt = read.table('./processed data/Final data/data_for_matching.txt')
dt = dt[!is.na(dt$avg_social_distance_score),]
dt = dt[complete.cases(dt),]
dt$bl_hisp = dt$hispanic + dt$black

attach(dt)

dt$IV = dt$avg_social_distance_score

n_t = sum(dt$treatment)
n_c = dim(dt)[1] - n_t
X = cbind(female, above_65, median_house_income, unemployment, bl_hisp,
          smoking, flu_vaccination, some_college, social_association,
          traffic, black, hispanic, RUCC, poverty, population.density,
          avg_obs_ill_0316_0322, avg_exp_ill_0316_0322, avg_exp_ill_0330_0405)
detach(dt)

# Perform matching
select1row = seq(1, nrow(X), 1)
mindist = 0.7
sinks = 0.8
matching = nbp_match_reverse_caliper(select1row, X, dt$IV, mindist, sinks)

for (i in 1:dim(matching)[1]){

}


ind = which((as.numeric(as.character(matching$Group1.ID)) <= 3020) &
              (as.numeric(as.character(matching$Group2.ID)) <= 3020))
mc = matching[ind,]

enc_ind = NULL
control_ind = NULL
for (i in 1:dim(mc)[1]){
  #cat(i, '\n')
  id_1 = as.numeric(as.character(mc$Group1.ID[i]))
  id_2 = as.numeric(as.character(mc$Group2.ID[i]))
  if (dt$avg_social_distance_score[id_1] > dt$avg_social_distance_score[id_2]){
    if (! id_1 %in% enc_ind) {
      enc_ind = c(enc_ind, id_1)
      control_ind = c(control_ind, id_2)
    }
  }
  if (dt$avg_social_distance_score[id_1] < dt$avg_social_distance_score[id_2]){
    if (! id_2 %in% enc_ind) {
        enc_ind = c(enc_ind, id_2)
        control_ind = c(control_ind, id_1)
  }
  }
}


dt_enc = dt[enc_ind,]
#dt_enc = dt_enc[!duplicated(dt_enc),]
dt_control = dt[control_ind,]
#dt_control = dt_control[!duplicated(dt_control),]
dt_matched = rbind(cbind(dt_enc, encouraged = rep(1, dim(dt_enc)[1])),
                   cbind(dt_control, encouraged = rep(0, dim(dt_control)[1])))

tb = xBalance(encouraged~female + above_65 + median_house_income+unemployment+
           smoking + flu_vaccination + some_college + social_association +
           traffic + bl_hisp + RUCC + poverty + population.density +
           avg_obs_ill_0316_0322 + avg_exp_ill_0316_0322 + avg_exp_ill_0330_0405 +
           avg_social_distance_score,
           data = dt_matched, report = c('adj.means', 'std.diffs', 'p.values'))

dim(dt_enc)[1]

plot_dt = data.frame(sd_score = dt_enc$avg_social_distance_score - dt_control$avg_social_distance_score)
ggplot(data = plot_dt, aes(x = sd_score)) + geom_histogram(binwidth = 0.1) + 
  xlab('matched pair difference in average social distancing score')


