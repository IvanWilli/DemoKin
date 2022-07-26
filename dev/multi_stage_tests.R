library(R.matlab)

# Program: matrix_construction_4867.m
# This Matlab script reads in tables giving life tables and parity transitions for Slovakia (obtained from the Human Mortality and Human Fertility databases, respectively), and constructs the age-specific parity transition matrices, the parity-specific survival matrices, the age-specific fertility matrices, and the parity-specific offspring assignment matrices.

matrices <- readMat("../../multi_stage_caswell/SVK_kinmats/SVKmats1952.mat")
names(matrices)

# omega
matrices$om

# stages
matrices$s

# stage transition prob: list of w elements with s x s matrix
matrices$U

# age advance (depends if U includes survival or not): list of s elements with w x w matrix
matrices$D

# state-specific fertility: list of w elements with s x s matrix
matrices$F

# assigns the offspring of individuals in stage j to the appropriate age class: list of s elements with w x w matrix
matrices$H

# Program: calling_kinship_SVK_4867.m
# This Matlab script reads in set of matrices produced by the program matrix_construction_4867.m and then calls the function kinship_function_parity_4867.m to obtain the resulting kinship network. The results are saved in the folder SVK_kinout with filenames SVKkinout1960.mat, etc.

years=1960:2014;

numyears=length(years);  #specific to SVK data
#add path to location of matrices
addpath('SVK_kinmats/')

for (iy in 1:numyears){

  iy = 1
  year=years[iy]

  #specify name of matrix file
  matrixs <- readMat(paste0("../../multi_stage_caswell/SVK_kinmats/SVKmats",year,".mat"))
  om = matrices$om
  s = matrices$s
  U = matrices$U
  D = matrices$D
  F = matrices$F
  H  = matrices$H

  #create the block diagonal matrices

  #identity matrices that are useful
  Iom=diag(1,om, om);
  Is=diag(1,s, s);

  bbU=matrix(0,s*om,s*om);
  bbF=matrix(0,s*om,s*om);
  for (i in 1:om){
    bbU = bbU + kron(Iom(:,i)*Iom(i,:),U{i}); #D is not useful
    bbF = bbF + kron(Iom(:,i)*Iom(i,:),F{i});
  }
  bbD=zeros(s*om);
  bbH=zeros(s*om);
  for (i in 1:s){
    bbD = bbD+kron(Is(:,i)*Is(i,:),D{i});
    bbH = bbH+kron(Is(:,i)*Is(i,:),H{i});
  }

  # create the age-stage matrices using the vec permuation formula
  K=vecperm(s,om);
  Ut= t(K)*bbD*K*bbU;
  Ft= t(K)*bbH*K*bbF;

  #conditional transition matrix, conditional on survival
  Gt=Ut*pinv(diag(sum(Ut)));

  #calculate distributions of mothers
  #projection matrix Atilde
  At=Ut+Ft;
  #eigenvalues and right eigenvectors
  [wt,d]=eig(At);
  d=diag(d);
  #find maximum eigenvalue
  pick=find(d==max(d));
  wt=wt(:,pick);
  #stable age-parity distribution normalized to sum to 1
  wt=wt/sum(wt);
  lambda=d(pick)

  #age-stage distribution of mothers
  pit=t(Ft(1,:)).*wt;
    pit=pit/sum(pit);
    # marginal age distribution of mothers
    piage=kron(Iom,t(ones(s,1))*pit;

#add path to call the kinship program
path('../',path)

#call the kinship function
kinout=kinship_function_parity(Ut,Ft,Gt,wt,pit,piage);

#save the kin output
#include path to output folder
myname=char(['SVK_kinout/SVKkinout' num2str(years(iy)) '.mat'])
save(myname,'kinout')

}


# Program: kinship_function_parity_4867.m
# This Matlab function takes as input a set of block-structured transition and fertility matrices and returns the kinship results.




