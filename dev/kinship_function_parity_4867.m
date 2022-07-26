function out=kinship_function_parity(Ut,Ft,Gt,wt,pit,piage)
%
%function to compute kinship network for multistate age x parity model
% Supplement to:
% Caswell, H. 2020. The formal demography of kinship II. Multistate models,
% parity, and sibship. Demographic Research 42:1097-1144
%
% Has been successfully used under Matlab R2018b
%
%
%inputs
% Ut=age-stage transition matrix
% Ft = age-stage fertility matrix
% Gt=age-stage transition matrix conditional on survival
% wt=stable age-stage distribution, normalized to sum to 1
% pit=age-stage distribution of mothers
% piage = marginal age distribution of mothers


%number of age classes
om=length(piage);
%number of stages
s=length(pit)/om;

%identity matrices useful in calculations
Iom=eye(om);
Is=eye(s);
Isom=eye(s*om);

%frequently used zero vector for initial condition
zvec=zeros(s*om,1);

%frequently used om-1 limit for iterations
omz=om-1;

% the following code calculates age-stage distributions, 
% for each type of kin, for each age x of Focal,
% and stores these as columns of an array
% e.g., a(x) = daughters at age x; A(:,x) contains a(x)

% dynamics of Focal
% initial condition
phiz=zeros(s*om,1);
phiz(1)=1;
%age-stage vector of Focal, conditional on survival
Phi(:,1)=phiz;
for ix=1:omz
    Phi(:,ix+1)=Gt*Phi(:,ix);
end

% a: daughters of focal

az=zvec;
A(:,1)=az;
for ix=1:omz
    A(:,ix+1)=Ut*A(:,ix) + Ft*Phi(:,ix);
end % for ix


% b = granddaughters of Focal
b=zvec;
B(:,1)=b;
for ix=1:omz
    B(:,ix+1)=Ut*B(:,ix) + Ft*A(:,ix);
end


% c = greatgranddaughters of Focal
c=zvec;
C(:,1)=c;
for ix=1:omz
    C(:,ix+1)=Ut*C(:,ix) +Ft*B(:,ix);
end


% d = mothers of Focal
% conditional on mother having parity >0 

%momarray is an array with pit in each column
momarray=pit*ones(1,om);

Z=eye(s);
Z(1,1)=0;
for imom=1:om  %go through all columns of momarray
    E=Iom(:,imom)*Iom(imom,:);
    momarray(:,imom)=kron(E,Z)*momarray(:,imom);
    %selects age imom, and eliminates the zero parity row of momarray
        
end
%rescale columns of momarray to sum to 1
momarray=momarray*pinv(diag(sum(momarray)));

%set dzero to the average of the momarray over the ages of moms at birth of
%children
dzero=momarray*piage;

D(:,1)=dzero;
for ix=1:omz
    D(:,ix+1)=Ut*D(:,ix);
end


% g = maternal grandmothers of Focal
gzero=D*piage;

G(:,1)=gzero;
for ix=1:omz
    G(:,ix+1)=Ut*G(:,ix);
end


% h = great-grandmothers of Focal
hzero=G*piage;
H(:,1)=hzero;
for ix=1:omz
    H(:,ix+1)=Ut*H(:,ix) + 0;
end

% m = older sisters of Focal
mzero=A*piage;
M(:,1)=mzero;
for ix=1:omz
    M(:,ix+1)=Ut*M(:,ix) + 0;
end

% n = younger sisters of Focal
nzero=zvec;
N(:,1)=nzero;
for ix=1:omz
    N(:,ix+1)=Ut*N(:,ix) + Ft*D(:,ix);
end


% p = nieces through older sisters of Focal
pzero=B*piage;
P(:,1)=pzero;
for ix=1:omz
    P(:,ix+1)=Ut*P(:,ix) + Ft*M(:,ix);
end

% q = nieces through younger sisters of Focal
qzero=zvec;
Q(:,1)=qzero;
for ix=1:omz
    Q(:,ix+1)=Ut*Q(:,ix) + Ft*N(:,ix);
end

% r = aunts older than mother of Focal
rzero=M*piage;
R(:,1)=rzero;
for ix=1:omz
    R(:,ix+1)=Ut*R(:,ix) + 0;
end

% s = aunts younger than mother of Focal
szero=N*piage;
S(:,1)=szero;
for ix=1:omz
    S(:,ix+1)=Ut*S(:,ix) + Ft*G(:,ix);
end

% t = cousins from aunts older than mother of Focal
tzero=P*piage;
T(:,1)=tzero;
for ix=1:omz
    T(:,ix+1)=Ut*T(:,ix) + Ft*R(:,ix);
end


% v = cousins from aunts younger than mother of Focal
vzero=Q*piage;
V(:,1)=vzero;
for ix=1:omz
    V(:,ix+1)=Ut*V(:,ix) + Ft*S(:,ix);
end %for i


%overall kinship matrices, concatenating all kin
allkin=cat(3,A,B,C,D,G,H,M,N,P,Q,R,S,T,V);

%combining older and younger categories
% for sisters, neices, aunts, and cousins
allkin2=cat(3,A,B,C,D,G,H,M+N,P+Q,R+S,T+V);

%output structure
out.allkin=allkin;
out.allkin2=allkin2;
out.Phi=Phi;
out.pit=pit;
out.piage=piage;
out.om=om;
out.s=s;
 out.Ut=Ut;
out.Ft=Ft;
out.Gt=Gt;




