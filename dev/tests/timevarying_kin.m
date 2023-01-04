function kout=timevarying_kin(U,F,pi,om,pkin);
% function to return kinship network
% calculated from the rates and the kinship at the previous time
% U=survival matrix
% F=fertility matrix
% pi = distribution of ages of mothers
% om=number of age classes
% pkin = the array of all kin from the previous time step
% model structure
% k(x+1,t+1)=U(t)*k(x,t) + F(t)*kstar(x,t) for some other kin kstar

%set to full in case they arrive as sparse matrices
U=full(U);
F=full(F);
pi=full(pi);

%frequently used zero vector for initial condition
zvec=zeros(om,1);
I=eye(om);
omz=om-1;

% a: daughters of focal

A(:,1)=zvec;
for ix=1:omz
    ap=U*pkin.A(:,ix) + F*I(:,ix);
    A(:,ix+1)=ap;
    
end % for ix

% b = granddaughters of Focal

B(:,1)=zvec;
for ix=1:omz
    bp=U*pkin.B(:,ix) + F*pkin.A(:,ix);
    B(:,ix+1)=bp;
    
end


% c = greatgranddaughters of Focal
C(:,1)=zvec;
for ix=1:omz
    cp=U*pkin.C(:,ix) +F*pkin.B(:,ix);
    C(:,ix+1)=cp;
    
end


% d = mothers of Focal
D(:,1)=pi;
for ix=1:omz
    dp=U*pkin.D(:,ix) + 0;
    D(:,ix+1)=dp;
    
end


% g = grandmothers of Focal
%only maternal grandmothers right now
G(:,1)=pkin.D*pi;;
for ix=1:omz
    gp=U*pkin.G(:,ix) + 0;
    G(:,ix+1)=gp;
    
end


% h = greattrandmothers of Focal

H(:,1)=pkin.G*pi;
for ix=1:omz
    hp=U*pkin.H(:,ix) + 0;
    H(:,ix+1)=hp;
    
end

% m = older sisters of Focal

M(:,1)=pkin.A*pi;
for ix=1:omz
    mp=U*pkin.M(:,ix) + 0;
    M(:,ix+1)=mp;
    
end

% n = younger sisters

N(:,1)=zvec;
for ix=1:omz
    np=U*pkin.N(:,ix) + F*pkin.D(:,ix);
    N(:,ix+1)=np;
    
end


% p = nieces through older sisters

P(:,1)=pkin.B*pi;
for ix=1:omz
    pp=U*pkin.P(:,ix) + F*pkin.M(:,ix);
    P(:,ix+1)=pp;
end

% q = nieces through younger sisters

Q(:,1)=zvec;
for ix=1:omz
    qp=U*pkin.Q(:,ix) + F*pkin.N(:,ix);
    Q(:,ix+1)=qp;
    
end

% r = aunts older than mother

R(:,1)=pkin.M*pi;
for ix=1:omz
    rp=U*pkin.R(:,ix) + 0;
    R(:,ix+1)=rp;
    
end

% s = aunts younger than mother

S(:,1)=pkin.N*pi;
for ix=1:omz
    sp=U*pkin.S(:,ix) + F*pkin.G(:,ix);
    S(:,ix+1)=sp;
    
end

% t = cousins from older aunts

T(:,1)=pkin.P*pi;
for ix=1:omz
    tp=U*pkin.T(:,ix) + F*pkin.R(:,ix);
    T(:,ix+1)=tp;
    
end


% v = cousins from aunts younger than mother

V(:,1)=pkin.Q*pi;
for ix=1:omz
    vp=U*pkin.V(:,ix) + F*pkin.S(:,ix);
    V(:,ix+1)=vp;
    
end

%concatenate kin matrices
allkin=cat(3,A,B,C,D,G,H,M,N,P,Q,R,S,T,V);

%concatenate, combining older and younger sisters, etc.
allkin2=cat(3,A,B,C,D,G,H,M+N,P+Q,R+S,T+V);

%create output structures
kout.A=A;
kout.B=B;
kout.C=C;
kout.D=D;
kout.G=G;
kout.H=H;
kout.M=M;
kout.N=N;
kout.P=P;
kout.Q=Q;
kout.R=R;
kout.S=S;
kout.T=T;
kout.V=V;

kout.allkin=allkin;
kout.allkin2=allkin2;

kout.U=U;
kout.F=F;
kout.pi=pi;

 