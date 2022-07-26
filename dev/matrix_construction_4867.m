
% script to prepare matrices for multistate age x parity model
% Supplement to:
% Caswell, H. 2020. The formal demography of kinship II. Multistate models,
% parity, and sibship. Demographic Research 42:1097-1144
%
%requires Matlab Table files obtained from HMD (fltper) and HFD (mi)
% Has been successfully used under Matlab R2018b

%add folder contraining the table files to path
addpath('SVK_tables/')

% load the female lifetable file
load('SVKfltperTable.mat')
%columns in this Table: Year,Age,mx,qx,ax,lx,dx,Lx,Tx,ex
lt=ltable;

%load the parity state transition file
load('SVKmiTable.mat')
%columns: Year,Age,mi1,mi2,mi3,mi4,mi5p

%find year ranges
minfertyear=min(fert.Year);
maxfertyear=max(fert.Year);

minltyear=min(lt.Year);
maxltyear=max(lt.Year);

%pick a starting year and ending year
startyear=max([minfertyear minltyear]);
endyear=min([maxfertyear maxltyear]);

%array of years and number of years
years=startyear:endyear;
numyears=endyear-startyear+1;

for iy=1:numyears
    years(iy);
    
    %find life table and qx array for year iy
    pick=find(lt.Year==years(iy));
    qx=table2array(lt(pick,4));
    
    %find fertility and create fertility array
    pick=find(fert.Year==years(iy));
    fertarray=table2array(fert(pick,[2:7]));
    
    %number of age classes
    %om=length(qx)-1;
    om=length(qx)-1;
    %number of parity classes
    s=6;
    
    %extend the fertility array
    startfert=fertarray(1,1);
    endfert=fertarray(end,1);
    %put zeros before age of first reproduction
    fertarray=[zeros(startfert-1,6); fertarray];
    fertarray(1:startfert-1,1)=(1:startfert-1)';
    %put zeros after age of last reproduction
    fertarray=[fertarray; zeros(om-endfert,6)];
    fertarray(endfert+1:om,1)=(endfert+1:om)';
    
    %remove age column from fertarray
    fertarray=fertarray(:,2:6);
    
    %construct the stage transition matrices using probabilities
    for i=1:om
        U{i} = diag(fertarray(i,:),-1);
        %transform subdiagonals to probabilities
        U{i}=U{i}./(1+0.5*U{i});
        %fill in diagonal entries
        U{i}=U{i}+diag([1-diag(U{i},-1) ; 1]);
    end
    
    %construct the age transition and survival matrices
    for i=1:s
        D{i}=diag(1-qx(1:om-1),-1);
    end
    
    %construct fertility matrices
    for i=1:om
        F{i}=zeros(s,s);
        F{i}(1,1:s-1)=diag(U{i},-1);
        F{i}(1,s)=U{i}(s,s-1);
        %divide fertility by 2
        F{i}=F{i}/2;
    end
    
    %stage assignment matrices
    for i=1:s
        H{i}=zeros(om,om);
        H{i}(1,:)=1;
    end
    
    %include path to folder where matrix files are to be stored
    myname=char(['SVK_kinmats/SVKmats' num2str(years(iy)) '.mat'])
    %save the matrices into a .mat file
    save(myname,'U','D','F','H','om','s')
    
end

