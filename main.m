% This is the main code to solve and simluate the Airbnb model.
clear all

%% Algorithm parmeters
n      = 5; %number of control 
ne     = 4; %number of exogenous shocks
npoly  = 4; %degree of polynomial
ncoef  = 1+4*n+2*n*(n-1); % # of polynomials in Smolyak
nnode  = npoly+1; % # of Nodes
nstate = 7; %number of Gauss-Hermite divisions.

%% Structural parameters
alpha  = 0.8;
beta   = 0.9;
sigma  = 5;
meanRS = 1.08;
meanRB = 1.02;
meanph = 1;
meanpe = 0.0439;
meanpa = 0.0698;
rhoy   = 0.748;
rhoph  = 0.7960;
rhope  = 0.8341;
rhopa  = 0.7503;
ses    = 0.158;
sey    = 0.15;
seph   = 0.0424;
sepe   = 0.0510;
sepa   = 0.1559;

%% Creating Chebyshev extrema nodes
set1    = chebyextrema(1);
set2    = chebyextrema(2);
set3    = chebyextrema(3);
set4    = chebyextrema(4);
set2(2) = 0;
set3(3) = 0;
set4(5) = 0;
% reorder to Smolyak
set2([1 2]) = set2([2 1]); 
set3([1 2 3 4 5]) = set3([3 1 5 2 4 ]);

%% Build sparse grid
G11111 = cartesianproduct5(set1,set1,set1,set1,set1);
G11112 = cartesianproduct5(set1,set1,set1,set1,set2);
G11113 = cartesianproduct5(set1,set1,set1,set1,set3);
G11121 = cartesianproduct5(set1,set1,set1,set2,set1);
G11122 = cartesianproduct5(set1,set1,set1,set2,set2);
G11131 = cartesianproduct5(set1,set1,set1,set3,set1);
G11211 = cartesianproduct5(set1,set1,set2,set1,set1);
G11212 = cartesianproduct5(set1,set1,set2,set1,set2);
G11221 = cartesianproduct5(set1,set1,set2,set2,set1);
G11311 = cartesianproduct5(set1,set1,set3,set1,set1);
G12111 = cartesianproduct5(set1,set2,set1,set1,set1);
G12112 = cartesianproduct5(set1,set2,set1,set1,set2);
G12121 = cartesianproduct5(set1,set2,set1,set2,set1);
G12211 = cartesianproduct5(set1,set2,set2,set1,set1);
G13111 = cartesianproduct5(set1,set3,set1,set1,set1);
G21111 = cartesianproduct5(set2,set1,set1,set1,set1);
G21112 = cartesianproduct5(set2,set1,set1,set1,set2);
G21121 = cartesianproduct5(set2,set1,set1,set2,set1);
G21211 = cartesianproduct5(set2,set1,set2,set1,set1);
G22111 = cartesianproduct5(set2,set2,set1,set1,set1);
G31111 = cartesianproduct5(set3,set1,set1,set1,set1);
G   = [G11111;G11112;G11113;G11121;G11122;G11131;G11211;G11212;G11221;G11311;G12111;G12112;G12121;G12211;G13111;G21111;G21112;G21121;G21211;G22111;G31111];
G75 = unique(G,'row');
index=zeros(ncoef,n);
%find index for each node.
for i=1:ncoef
    for j=1:n
        index(i,j) = find(set3==G75(i,j));
    end
end
indexth = [1 1 1 1 1;
    1 1 1 1 2;
    1 1 1 1 3;
    1 1 1 2 1;
    1 1 1 3 1;
    1 1 2 1 1;
    1 1 3 1 1;
    1 2 1 1 1;
    1 3 1 1 1;
    2 1 1 1 1;
    3 1 1 1 1;
    1 1 1 1 4;
    1 1 1 1 5;
    1 1 1 2 2;
    1 1 1 2 3;
    1 1 1 3 2;
    1 1 1 3 3;
    1 1 1 4 1;
    1 1 1 5 1;
    1 1 2 1 2;
    1 1 2 1 3;
    1 1 3 1 2;
    1 1 3 1 3;
    1 1 2 2 1;
    1 1 2 3 1;
    1 1 3 2 1;
    1 1 3 3 1;
    1 1 4 1 1;
    1 1 5 1 1;
    1 2 1 1 2;
    1 2 1 1 3;
    1 3 1 1 2;
    1 3 1 1 3;
    1 2 1 2 1;
    1 2 1 3 1;
    1 3 1 2 1;
    1 3 1 3 1;
    1 2 2 1 1;
    1 2 3 1 1;
    1 3 2 1 1;
    1 3 3 1 1;
    1 4 1 1 1;
    1 5 1 1 1;
    2 1 1 1 2;
    2 1 1 1 3;
    3 1 1 1 2;
    3 1 1 1 3;
    2 1 1 2 1;
    2 1 1 3 1;
    3 1 1 2 1;
    3 1 1 3 1;
    2 1 2 1 1 ;
    2 1 3 1 1 ;
    3 1 2 1 1;
    3 1 3 1 1;
    2 2 1 1 1;
    2 3 1 1 1;
    3 2 1 1 1;
    3 3 1 1 1;
    4 1 1 1 1;
    5 1 1 1 1];

%% Initialize hermite nodes for approximation:log(y), log(ph), and log(pe)
[wgridy,wmaty]    = hernodes(nstate);
wgridy            = wgridy*sqrt(2)*sey;
[wgridph,wmatph]  = hernodes(nstate);
wgridph           = wgridph*sqrt(2)*seph;
[wgridpe,wmatpe]  = hernodes(nstate);
wgridpe           = wgridpe*sqrt(2)*sepe;
[wgridpa,wmatpa]  = hernodes(nstate);
wgridpa           = wgridpa*sqrt(2)*sepa;
weight = kron(kron(kron(wmatph,wmaty),wmatpe),wmatpa);

%% Transform the domain of state variables and create cheby poly. 
%log(fw), not fw
fwmin    = log(0.4);
fwmax    = log(4);
fwt      = itransfo(set3,fwmin,fwmax); %[fwmin,fwmax]
MFW      = cheb(set3',0:npoly);%poly based on extrema,row is number of nodes
%log(y), not y.
ymin     = wgridy(nstate);
ymax     = wgridy(1);
yt       = itransfo(set3,ymin,ymax);%[-1,1] to [ymin, ymax]
MY       = cheb(set3',0:npoly);
%log(ph), not ph.
phmin    = log(meanph)-3*seph;
phmax    = log(meanph)+3*seph;
pht      = itransfo(set3,phmin,phmax);
MPH      = cheb(set3',0:npoly); %
%log(pe)
pemin    = log(meanpe)-3*sepe;
pemax    = log(meanpe)+3*sepe;
pet      = itransfo(set3,pemin,pemax);
MPE      = cheb(set3',0:npoly); %
pamin    = log(meanpa)-3*sepa;
pamax    = log(meanpa)+3*sepa;
pat      = itransfo(set3,pamin,pamax);
MPA      = cheb(set3',0:npoly); %
M        = MY;

%% Initial condition
a0 = [zeros(2*ncoef,1);zeros(ncoef,1);0.54;zeros(ncoef-1,1);0.03;zeros(ncoef-1,1)];
a0 = a0(:);

%% Main loop
tic
aparam    = [n,ne,npoly,ncoef,nnode,nstate]';
sparam    = [alpha beta sigma meanRS meanRB rhoy rhoph rhope rhopa meanph meanpe meanpa]';
options   = optimoptions('lsqnonlin','Display','iter','MaxFunctionEvaluations',300,'MaxIterations',1000,'FunctionTolerance',1e-6,'UseParallel',true);
LB        = [repelem(-Inf,ncoef*2),repelem(-10,ncoef*3)];
UB        = [repelem(Inf,ncoef*2),repelem(10,ncoef*3)];
[th,fval] = lsqnonlin('residuals',a0,LB,UB,options,aparam,sparam,fwt,yt,pht,pet,pat,wgridy,wgridpe,wgridph,wgridpa,weight,index,indexth,M); %step 6
toc

%% Simulation
T = simulation(th,aparam,sparam,indexth,sey,seph,sepe,sepa,fwmin,fwmax,phmin,phmax,ymin,ymax,pemin,pemax,pamin,pamax,M);
T


