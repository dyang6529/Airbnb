% Simulation

function T = simulation(th,aparam,sparam,indexth,sey,seph,sepe,sepa,fwmin,fwmax,phmin,phmax,ymin,ymax,pemin,pemax,pamin,pamax,M)
n      = aparam(1);
nex    = aparam(2);
npoly  = aparam(3);
ncoef  = aparam(4);
nnode  = aparam(5);
nstate = aparam(6);
alpha  = sparam(1);
beta   = sparam(2);
sigma  = sparam(3);
meanRS = sparam(4);
meanRB = sparam(5);
rhoy   = sparam(6);
rhoph  = sparam(7);
rhope  = sparam(8);
rhopa  = sparam(9);
meanph = sparam(10);
meanpe = sparam(11);
meanpa = sparam(12);

policy_c  = th(1:ncoef);
policy_h  = th(ncoef+1:2*ncoef);
policy_b  =  th(2*ncoef+1:3*ncoef);
policy_dh = th(3*ncoef+1:4*ncoef);
policy_da = th(4*ncoef+1:end);

%% Initialization
T           = 100000;
consumption = zeros(T,1)';
bond        = zeros(T,1)';
house       = zeros(T,1)';
deltah      = zeros(T,1)';
deltaa      = zeros(T,1)';
fin_wealth  = zeros(T,1)';
y           = zeros(T,1)';
ph          = zeros(T,1)';
pe          = zeros(T,1)';
pa          = zeros(T,1)';
% Inital values
bond(1)     = 0;
house(1)    = 1;
deltah(1)   = 0;
deltaa(1)   = 0;
y(1)        = 0;
ph(1)       = 0;
pe(1)       = -3;
pa(1)       = -2.5257;
%% Simulation
tic
for t = 1:T
    fin_wealth(t) = log(meanRB*bond(t)+(1-deltah(t)-deltaa(t))*exp(pe(t))*house(t)+exp(ph(t))*house(t)+deltaa(t)*exp(pa(t))*house(t)); 
    root_fwt1     = transfo(fin_wealth(t),min(fwmin,fin_wealth(t)),max(fin_wealth(t),fwmax));
    MFW1          = cheb(root_fwt1,0:npoly);
    root_yt1      = transfo(y(t),min(ymin,y(t)),max(ymax,y(t)));
    MY1           = cheb(root_yt1,0:npoly);
    root_pht1     = transfo(ph(t),min(phmin,ph(t)),max(phmax,ph(t))); 
    MPH1          = cheb(root_pht1,0:npoly);
    root_pet1     = transfo(pe(t),min(pemin,pe(t)),max(pemax,pe(t))); 
    MPE1          = cheb(root_pet1,0:npoly);
    root_pat1     = transfo(pa(t),min(pamin,pa(t)),max(pamax,pa(t))); 
    MPA1          =  cheb(root_pat1,0:npoly);
    V0=zeros(ncoef,1)';
    for k = 1:ncoef
        V0(k) = MFW1(indexth(k,1))*MY1(indexth(k,2))*MPH1(indexth(k,3))*MPE1(indexth(k,4))*MPA1(indexth(k,5));
    end
    consumption(t) = exp(V0*policy_c);
    house(t+1)     = exp(V0*policy_h);
    bond(t+1)      = max(V0*policy_b,0)^2;
    deltah(t+1)    = min(max(V0*policy_dh,0)^2,1);
    deltaa(t+1)    = min(max(V0*policy_da,0)^2,1);
    y(t+1)         = rhoy*y(t)+normrnd(0,sey); 
    ph(t+1)        = rhoph*ph(t)+(1-rhoph)*log(meanph)+normrnd(0,seph);
    pe(t+1)        = rhope*pe(t)+(1-rhope)*log(meanpe)+normrnd(0,sepe);
    pa(t+1)        = rhopa*pa(t)+(1-rhopa)*log(meanpa)+normrnd(0,sepa);
end
toc
% Burner
consumption    =consumption(5000:end);
house          =house(5000:end);
bond           =bond(5000:end);
fin_wealth     =fin_wealth(5000:end);
y              =y(5000:end-1);
ph             =ph(5000:end-1);
pe             =pe(5000:end-1);
pa             =pa(5000:end-1);
deltah         =deltah(5000:end);
deltaa         =deltaa(5000:end);

%% Find moments
all        = [consumption;house(2:end);bond(2:end);fin_wealth;y;ph;deltah(2:end);pe;deltaa(2:end);pa;house(2:end).*deltah(2:end);house(2:end).*deltaa(2:end)]';
Mcorr      = corrcoef(all)
house_avg  = mean(house(1:end-1).*exp(ph))/mean(exp(fin_wealth)-exp(pe).*(1-deltah(1:end-1)-deltaa(1:end-1)).*house(1:end-1)-deltaa(1:end-1).*exp(pa).*house(1:end-1));
bond_avg   = mean(meanRB*bond(1:end-1))/mean(exp(fin_wealth)-exp(pe).*(1-deltah(1:end-1)-deltaa(1:end-1)).*house(1:end-1)-deltaa(1:end-1).*exp(pa).*house(1:end-1));
deltah_avg = mean(deltah);
deltaa_avg = mean(deltaa);
phminimum  = exp(min(ph));
phmaximum  = exp(max(ph));
house_size = mean(house);
total_weal = exp(mean(fin_wealth));
avg_ph     = exp(mean(ph));
Rental_wea = mean(exp(pe).*(1-deltah(1:end-1)-deltaa(1:end-1)).*house(1:end-1))/mean(exp(fin_wealth));
Airbnb_wea = mean(deltaa(1:end-1).*exp(pa).*house(1:end-1))/mean(exp(fin_wealth));
consumpt   = mean(consumption);
house_cons = mean(house.*deltah);
house_var  = var(house./house_avg);
deltah_var = var(deltah./deltah_avg);
houseadj   = sum(house(1:end-1)==house(2:end))/length(house);
% Checking boundary violations
fwUB       = sum(fin_wealth>fwmax)/length(fin_wealth); %less than 1% violated
fwLB       = sum(fin_wealth<fwmin)/length(fin_wealth); %less than 1% violated
YUB        = sum(y>ymax)/length(fin_wealth); %none violated
YLB        = sum(y<ymin)/length(fin_wealth); %none violated
% Table
Name  = {'Avg Housing Wealth';'Average Bond Wealth';'Delta H';'Delta A';'FW upperbound';'FW lowerbound';'Y UB';'Y LB';'Ph minimum';'Ph maximum';'Average ph';...
    'House size';'Financial wealth';'Rental wealth';'Airbnb_wealth';'Numeraire Consumption';'Housing Consumption';'Normalized Variance for Housing';...
    'Normalized Variance for DeltaH';'Same house proportion'};
Value = [house_avg;bond_avg;deltah_avg;deltaa_avg;fwUB;fwLB;YUB;YLB;phminimum;phmaximum;avg_ph;house_size;total_weal;Rental_wea;Airbnb_wea;...
    consumpt;house_cons;house_var;deltah_var;houseadj];
T     = table(Name,Value);
end