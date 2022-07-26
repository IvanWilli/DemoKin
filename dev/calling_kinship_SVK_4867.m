%script to calculate kinship results
%this script calls the function kinship_function_parity_4867
%requires the function vecperm.m to create vec-permutation matrix
%
% Supplement to:
% Caswell, H. 2020. The formal demography of kinship II. Multistate models,
% parity, and sibship. Demographic Research 42:1097-1144
%
% Has been successfully used under Matlab R2018b

%specify range of years to analyze
years=1960:2014;

%years=2002;

numyears=length(years);  %specific to SVK data
%add path to location of matrices
addpath('SVK_kinmats/')

for iy=1:numyears
    year=years(iy)
    
    %specify name of matrix file
    fname=char(['SVKmats' num2str(1950+iy-1) '.mat']);
    %load matrix file
    load(fname)
    
    %create the block diagonal matrices
    
    %identity matrices that are useful
    Iom=eye(om);
    Is=eye(s);
    
    bbU=zeros(s*om);
    bbF=zeros(s*om);
    for i=1:om
        bbU = bbU + kron(Iom(:,i)*Iom(i,:),U{i});
        bbF = bbF + kron(Iom(:,i)*Iom(i,:),F{i});
    end
    bbD=zeros(s*om);
    bbH=zeros(s*om);
    for i=1:s
        bbD = bbD+kron(Is(:,i)*Is(i,:),D{i});
        bbH = bbH+kron(Is(:,i)*Is(i,:),H{i});
    end
    
    %create the age-stage matrices using the vec permuation formula
    K=vecperm(s,om);
    Ut= K'*bbD*K*bbU;
    Ft= K'*bbH*K*bbF;
    
    %conditional transition matrix, conditional on survival
    Gt=Ut*pinv(diag(sum(Ut)));
    
    %calculate distributions of mothers
    %projection matrix Atilde
    At=Ut+Ft;
    %eigenvalues and right eigenvectors
    [wt,d]=eig(At);
    d=diag(d);
    %find maximum eigenvalue
    pick=find(d==max(d));
    wt=wt(:,pick);
    %stable age-parity distribution normalized to sum to 1
    wt=wt/sum(wt);
    lambda=d(pick)
    
    %age-stage distribution of mothers
    pit=Ft(1,:)'.*wt;
    pit=pit/sum(pit);
    %marginal age distribution of mothers
    piage=kron(Iom,ones(s,1)')*pit;
    
    clear At
    
    %add path to call the kinship program
    path('../',path)
    
    %call the kinship function
    kinout=kinship_function_parity(Ut,Ft,Gt,wt,pit,piage);
    
    %save the kin output
    %include path to output folder
    myname=char(['SVK_kinout/SVKkinout' num2str(years(iy)) '.mat'])
    save(myname,'kinout')
end
