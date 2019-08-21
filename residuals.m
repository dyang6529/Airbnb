% Return a vector of residuals
function res = residuals(theta,aparam,sparam,fwt,yt,pht,pet,pat,wgridy,wgridpe,wgridph,wgridpa,weight,index,indexth,M)

n      = aparam(1);
ne    = aparam(2);
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

RHS1   = zeros(ncoef,1);
RHS2   = zeros(ncoef,1);
RHS3   = zeros(ncoef,1);
RHS4   = zeros(ncoef,1);
RHS5   = zeros(ncoef,1);
LHS1   = zeros(ncoef,1);
LHS2   = zeros(ncoef,1);
LHS3   = zeros(ncoef,1);
LHS4   = zeros(ncoef,1);
LHS5   = zeros(ncoef,1);

theta11 = theta(1:ncoef);
theta12 = theta(ncoef+1:2*ncoef);
theta13 = theta(2*ncoef+1:3*ncoef);
theta14 = theta(3*ncoef+1:4*ncoef);
theta15 = theta(4*ncoef+1:end);

%i for all nodes, j for all evaulations
dj  = zeros(ncoef,1);
daj = zeros(ncoef,1);
dbj = zeros(ncoef,1);
dcj = zeros(ncoef,1);
ddj = zeros(ncoef,1);
fj  = zeros(nnode,nnode,nnode,nnode,nnode);
faj = zeros(nnode,nnode,nnode,nnode,nnode);
fbj = zeros(nnode,nnode,nnode,nnode,nnode);
fcj = zeros(nnode,nnode,nnode,nnode,nnode);
fdj = zeros(nnode,nnode,nnode,nnode,nnode);

% Evaluate function d,da,db,dc,dd at all points, using inital guess theta
for j = 1:ncoef
    % Evalulate d using all theta
    djk  = zeros(ncoef,1);
    dajk = zeros(ncoef,1);
    dbjk = zeros(ncoef,1);
    dcjk = zeros(ncoef,1);
    ddjk = zeros(ncoef,1);
    for k = 1:ncoef
        djk(k)  = theta11(k)*M(index(j,1),indexth(k,1))*M(index(j,2),indexth(k,2))*M(index(j,3),indexth(k,3))*M(index(j,4),indexth(k,4))*M(index(j,5),indexth(k,5));
        dajk(k) = theta12(k)*M(index(j,1),indexth(k,1))*M(index(j,2),indexth(k,2))*M(index(j,3),indexth(k,3))*M(index(j,4),indexth(k,4))*M(index(j,5),indexth(k,5));
        dbjk(k) = theta13(k)*M(index(j,1),indexth(k,1))*M(index(j,2),indexth(k,2))*M(index(j,3),indexth(k,3))*M(index(j,4),indexth(k,4))*M(index(j,5),indexth(k,5));
        dcjk(k) = theta14(k)*M(index(j,1),indexth(k,1))*M(index(j,2),indexth(k,2))*M(index(j,3),indexth(k,3))*M(index(j,4),indexth(k,4))*M(index(j,5),indexth(k,5));
        ddjk(k) = theta15(k)*M(index(j,1),indexth(k,1))*M(index(j,2),indexth(k,2))*M(index(j,3),indexth(k,3))*M(index(j,4),indexth(k,4))*M(index(j,5),indexth(k,5));
    end
    dj(j)   = sum(djk);
    daj(j)  = sum(dajk);
    dbj(j)  = sum(dbjk);
    dcj(j)  = sum(dcjk);
    ddj(j)  = sum(ddjk);
    fj(index(j,1),index(j,2),index(j,3),index(j,4),index(j,5))  = dj(j);
    faj(index(j,1),index(j,2),index(j,3),index(j,4),index(j,5)) = daj(j);
    fbj(index(j,1),index(j,2),index(j,3),index(j,4),index(j,5)) = dbj(j);
    fcj(index(j,1),index(j,2),index(j,3),index(j,4),index(j,5)) = dcj(j);
    fdj(index(j,1),index(j,2),index(j,3),index(j,4),index(j,5)) = ddj(j);
end

%     Find the coefficients using function values
cl      = [1 2 2 1 1];
klin    = [1 3 3 5 5]; %inside summation
klout   = [2 3 3 5 5]-1; %outside summation
theta21 = zeros(ncoef,1)';
theta22 = zeros(ncoef,1)';
theta23 = zeros(ncoef,1)';
theta24 = zeros(ncoef,1)';
theta25 = zeros(ncoef,1)';
%theta 2 to 41.
for m = 2:ncoef
    cstth = 2^n/(cl(indexth(m,1))*cl(indexth(m,2))*cl(indexth(m,3))*cl(indexth(m,4))*cl(indexth(m,5))*klout(indexth(m,1))*klout(indexth(m,2))...
        *klout(indexth(m,3))*klout(indexth(m,4))*klout(indexth(m,5)));
    for i_th = 1:klin(indexth(m,1))
        for j_th = 1:klin(indexth(m,2))
            for k_th = 1:klin(indexth(m,3))
                for l_th = 1:klin(indexth(m,4))
                    for m_th = 1:klin(indexth(m,5))
                        theta21(m) = theta21(m)+cstth*M(i_th,indexth(m,1))*M(j_th,indexth(m,2))*M(k_th,indexth(m,3))*M(l_th,indexth(m,4))*M(m_th,indexth(m,5))...
                            *fj(i_th,j_th,k_th,l_th,m_th)/((cl(i_th)*cl(j_th))*cl(k_th)*cl(l_th)*cl(m_th));
                        theta22(m) = theta22(m)+cstth*M(i_th,indexth(m,1))*M(j_th,indexth(m,2))*M(k_th,indexth(m,3))*M(l_th,indexth(m,4))*M(m_th,indexth(m,5))...
                            *faj(i_th,j_th,k_th,l_th,m_th)/((cl(i_th)*cl(j_th))*cl(k_th)*cl(l_th)*cl(m_th));
                        theta23(m) = theta23(m)+cstth*M(i_th,indexth(m,1))*M(j_th,indexth(m,2))*M(k_th,indexth(m,3))*M(l_th,indexth(m,4))*M(m_th,indexth(m,5))...
                            *fbj(i_th,j_th,k_th,l_th,m_th)/((cl(i_th)*cl(j_th))*cl(k_th)*cl(l_th)*cl(m_th));
                        theta24(m) = theta24(m)+cstth*M(i_th,indexth(m,1))*M(j_th,indexth(m,2))*M(k_th,indexth(m,3))*M(l_th,indexth(m,4))*M(m_th,indexth(m,5))...
                            *fcj(i_th,j_th,k_th,l_th,m_th)/((cl(i_th)*cl(j_th))*cl(k_th)*cl(l_th)*cl(m_th));
                        theta25(m) = theta25(m)+cstth*M(i_th,indexth(m,1))*M(j_th,indexth(m,2))*M(k_th,indexth(m,3))*M(l_th,indexth(m,4))*M(m_th,indexth(m,5))...
                            *fdj(i_th,j_th,k_th,l_th,m_th)/((cl(i_th)*cl(j_th))*cl(k_th)*cl(l_th)*cl(m_th));
                    end
                end
            end
        end
    end
end
theta21(1)  = fj(1,1,1,1,1)-theta21(2)*M(1,2)-theta21(3)*M(1,3)-theta21(4)*M(1,2)-theta21(5)*M(1,3)-theta21(6)*M(1,2)-theta21(7)*M(1,3)-theta21(8)*M(1,2)-...
    theta21(9)*M(1,3)-theta21(10)*M(1,2)-theta21(11)*M(1,3)-theta21(12)*M(1,4)-theta21(13)*M(1,5)-theta21(14)*M(1,2)*M(1,2)-theta21(15)*M(1,2)*M(1,3)-...
    theta21(16)*M(1,3)*M(1,2)-theta21(17)*M(1,3)*M(1,3)-theta21(18)*M(1,4)-theta21(19)*M(1,5)-theta21(20)*M(1,2)*M(1,2)-theta21(21)*M(1,2)*M(1,3)-...
    theta21(22)*M(1,3)*M(1,2)-theta21(23)*M(1,3)*M(1,3)-theta21(24)*M(1,2)*M(1,2)-theta21(25)*M(1,2)*M(1,3)-theta21(26)*M(1,3)*M(1,2)-...
    theta21(27)*M(1,3)*M(1,3)-theta21(28)*M(1,4)-theta21(29)*M(1,5)-theta21(30)*M(1,2)*M(1,2)-theta21(31)*M(1,2)*M(1,3)-theta21(32)*M(1,3)*M(1,2)-...
    theta21(33)*M(1,3)*M(1,3)-theta21(34)*M(1,2)*M(1,2)-theta21(35)*M(1,2)*M(1,3)-theta21(36)*M(1,3)*M(1,2)-theta21(37)*M(1,3)*M(1,3)-...
    theta21(38)*M(1,2)*M(1,2)-theta21(39)*M(1,2)*M(1,3)-theta21(40)*M(1,3)*M(1,2)-theta21(41)*M(1,3)*M(1,3)-theta21(42)*M(1,4)-theta21(43)*M(1,5)-...
    theta21(44)*M(1,2)*M(1,2)-theta21(45)*M(1,2)*M(1,3)-theta21(46)*M(1,3)*M(1,2)-theta21(47)*M(1,3)*M(1,3)-theta21(48)*M(1,2)*M(1,2)-...
    theta21(49)*M(1,2)*M(1,3)-theta21(50)*M(1,3)*M(1,2)-theta21(51)*M(1,3)*M(1,3)-theta21(52)*M(1,2)*M(1,2)-theta21(53)*M(1,2)*M(1,3)-...
    theta21(54)*M(1,3)*M(1,2)-theta21(55)*M(1,3)*M(1,3)-theta21(56)*M(1,2)*M(1,2)-theta21(57)*M(1,2)*M(1,3)-theta21(58)*M(1,3)*M(1,2)-...
    theta21(59)*M(1,3)*M(1,3)-theta21(60)*M(1,4)-theta21(61)*M(1,5);
theta22(1)  = faj(1,1,1,1,1)-theta22(2)*M(1,2)-theta22(3)*M(1,3)-theta22(4)*M(1,2)-theta22(5)*M(1,3)-theta22(6)*M(1,2)-theta22(7)*M(1,3)-theta22(8)*M(1,2)-...
    theta22(9)*M(1,3)-theta22(10)*M(1,2)-theta22(11)*M(1,3)-theta22(12)*M(1,4)-theta22(13)*M(1,5)-theta22(14)*M(1,2)*M(1,2)-theta22(15)*M(1,2)*M(1,3)-...
    theta22(16)*M(1,3)*M(1,2)-theta22(17)*M(1,3)*M(1,3)-theta22(18)*M(1,4)-theta22(19)*M(1,5)-theta22(20)*M(1,2)*M(1,2)-theta22(21)*M(1,2)*M(1,3)-...
    theta22(22)*M(1,3)*M(1,2)-theta22(23)*M(1,3)*M(1,3)-theta22(24)*M(1,2)*M(1,2)-theta22(25)*M(1,2)*M(1,3)-theta22(26)*M(1,3)*M(1,2)-...
    theta22(27)*M(1,3)*M(1,3)-theta22(28)*M(1,4)-theta22(29)*M(1,5)-theta22(30)*M(1,2)*M(1,2)-theta22(31)*M(1,2)*M(1,3)-theta22(32)*M(1,3)*M(1,2)-...
    theta22(33)*M(1,3)*M(1,3)-theta22(34)*M(1,2)*M(1,2)-theta22(35)*M(1,2)*M(1,3)-theta22(36)*M(1,3)*M(1,2)-theta22(37)*M(1,3)*M(1,3)-...
    theta22(38)*M(1,2)*M(1,2)-theta22(39)*M(1,2)*M(1,3)-theta22(40)*M(1,3)*M(1,2)-theta22(41)*M(1,3)*M(1,3)-theta22(42)*M(1,4)-theta22(43)*M(1,5)-...
    theta22(44)*M(1,2)*M(1,2)-theta22(45)*M(1,2)*M(1,3)-theta22(46)*M(1,3)*M(1,2)-theta22(47)*M(1,3)*M(1,3)-theta22(48)*M(1,2)*M(1,2)-...
    theta22(49)*M(1,2)*M(1,3)-theta22(50)*M(1,3)*M(1,2)-theta22(51)*M(1,3)*M(1,3)-theta22(52)*M(1,2)*M(1,2)-theta22(53)*M(1,2)*M(1,3)-...
    theta22(54)*M(1,3)*M(1,2)-theta22(55)*M(1,3)*M(1,3)-theta22(56)*M(1,2)*M(1,2)-theta22(57)*M(1,2)*M(1,3)-theta22(58)*M(1,3)*M(1,2)-...
    theta22(59)*M(1,3)*M(1,3)-theta22(60)*M(1,4)-theta22(61)*M(1,5);
theta23(1)  = fbj(1,1,1,1,1)-theta23(2)*M(1,2)-theta23(3)*M(1,3)-theta23(4)*M(1,2)-theta23(5)*M(1,3)-theta23(6)*M(1,2)-theta23(7)*M(1,3)-theta23(8)*M(1,2)-...
    theta23(9)*M(1,3)-theta23(10)*M(1,2)-theta23(11)*M(1,3)-theta23(12)*M(1,4)-theta23(13)*M(1,5)-theta23(14)*M(1,2)*M(1,2)-theta23(15)*M(1,2)*M(1,3)-...
    theta23(16)*M(1,3)*M(1,2)-theta23(17)*M(1,3)*M(1,3)-theta23(18)*M(1,4)-theta23(19)*M(1,5)-theta23(20)*M(1,2)*M(1,2)-theta23(21)*M(1,2)*M(1,3)-...
    theta23(22)*M(1,3)*M(1,2)-theta23(23)*M(1,3)*M(1,3)-theta23(24)*M(1,2)*M(1,2)-theta23(25)*M(1,2)*M(1,3)-theta23(26)*M(1,3)*M(1,2)-...
    theta23(27)*M(1,3)*M(1,3)-theta23(28)*M(1,4)-theta23(29)*M(1,5)-theta23(30)*M(1,2)*M(1,2)-theta23(31)*M(1,2)*M(1,3)-theta23(32)*M(1,3)*M(1,2)-...
    theta23(33)*M(1,3)*M(1,3)-theta23(34)*M(1,2)*M(1,2)-theta23(35)*M(1,2)*M(1,3)-theta23(36)*M(1,3)*M(1,2)-theta23(37)*M(1,3)*M(1,3)-...
    theta23(38)*M(1,2)*M(1,2)-theta23(39)*M(1,2)*M(1,3)-theta23(40)*M(1,3)*M(1,2)-theta23(41)*M(1,3)*M(1,3)-theta23(42)*M(1,4)-theta23(43)*M(1,5)-...
    theta23(44)*M(1,2)*M(1,2)-theta23(45)*M(1,2)*M(1,3)-theta23(46)*M(1,3)*M(1,2)-theta23(47)*M(1,3)*M(1,3)-theta23(48)*M(1,2)*M(1,2)-...
    theta23(49)*M(1,2)*M(1,3)-theta23(50)*M(1,3)*M(1,2)-theta23(51)*M(1,3)*M(1,3)-theta23(52)*M(1,2)*M(1,2)-theta23(53)*M(1,2)*M(1,3)-...
    theta23(54)*M(1,3)*M(1,2)-theta23(55)*M(1,3)*M(1,3)-theta23(56)*M(1,2)*M(1,2)-theta23(57)*M(1,2)*M(1,3)-theta23(58)*M(1,3)*M(1,2)-...
    theta23(59)*M(1,3)*M(1,3)-theta23(60)*M(1,4)-theta23(61)*M(1,5);
theta24(1)  = fcj(1,1,1,1,1)-theta24(2)*M(1,2)-theta24(3)*M(1,3)-theta24(4)*M(1,2)-theta24(5)*M(1,3)-theta24(6)*M(1,2)-theta24(7)*M(1,3)-theta24(8)*M(1,2)-...
    theta24(9)*M(1,3)-theta24(10)*M(1,2)-theta24(11)*M(1,3)-theta24(12)*M(1,4)-theta24(13)*M(1,5)-theta24(14)*M(1,2)*M(1,2)-theta24(15)*M(1,2)*M(1,3)-...
    theta24(16)*M(1,3)*M(1,2)-theta24(17)*M(1,3)*M(1,3)-theta24(18)*M(1,4)-theta24(19)*M(1,5)-theta24(20)*M(1,2)*M(1,2)-theta24(21)*M(1,2)*M(1,3)-...
    theta24(22)*M(1,3)*M(1,2)-theta24(23)*M(1,3)*M(1,3)-theta24(24)*M(1,2)*M(1,2)-theta24(25)*M(1,2)*M(1,3)-theta24(26)*M(1,3)*M(1,2)-...
    theta24(27)*M(1,3)*M(1,3)-theta24(28)*M(1,4)-theta24(29)*M(1,5)-theta24(30)*M(1,2)*M(1,2)-theta24(31)*M(1,2)*M(1,3)-theta24(32)*M(1,3)*M(1,2)-...
    theta24(33)*M(1,3)*M(1,3)-theta24(34)*M(1,2)*M(1,2)-theta24(35)*M(1,2)*M(1,3)-theta24(36)*M(1,3)*M(1,2)-theta24(37)*M(1,3)*M(1,3)-...
    theta24(38)*M(1,2)*M(1,2)-theta24(39)*M(1,2)*M(1,3)-theta24(40)*M(1,3)*M(1,2)-theta24(41)*M(1,3)*M(1,3)-theta24(42)*M(1,4)-theta24(43)*M(1,5)-...
    theta24(44)*M(1,2)*M(1,2)-theta24(45)*M(1,2)*M(1,3)-theta24(46)*M(1,3)*M(1,2)-theta24(47)*M(1,3)*M(1,3)-theta24(48)*M(1,2)*M(1,2)-...
    theta24(49)*M(1,2)*M(1,3)-theta24(50)*M(1,3)*M(1,2)-theta24(51)*M(1,3)*M(1,3)-theta24(52)*M(1,2)*M(1,2)-theta24(53)*M(1,2)*M(1,3)-...
    theta24(54)*M(1,3)*M(1,2)-theta24(55)*M(1,3)*M(1,3)-theta24(56)*M(1,2)*M(1,2)-theta24(57)*M(1,2)*M(1,3)-theta24(58)*M(1,3)*M(1,2)-...
    theta24(59)*M(1,3)*M(1,3)-theta24(60)*M(1,4)-theta24(61)*M(1,5);
theta25(1)  = fdj(1,1,1,1,1)-theta25(2)*M(1,2)-theta25(3)*M(1,3)-theta25(4)*M(1,2)-theta25(5)*M(1,3)-theta25(6)*M(1,2)-theta25(7)*M(1,3)-theta25(8)*M(1,2)-...
    theta25(9)*M(1,3)-theta25(10)*M(1,2)-theta25(11)*M(1,3)-theta25(12)*M(1,4)-theta25(13)*M(1,5)-theta25(14)*M(1,2)*M(1,2)-theta25(15)*M(1,2)*M(1,3)-...
    theta25(16)*M(1,3)*M(1,2)-theta25(17)*M(1,3)*M(1,3)-theta25(18)*M(1,4)-theta25(19)*M(1,5)-theta25(20)*M(1,2)*M(1,2)-theta25(21)*M(1,2)*M(1,3)-...
    theta25(22)*M(1,3)*M(1,2)-theta25(23)*M(1,3)*M(1,3)-theta25(24)*M(1,2)*M(1,2)-theta25(25)*M(1,2)*M(1,3)-theta25(26)*M(1,3)*M(1,2)-...
    theta25(27)*M(1,3)*M(1,3)-theta25(28)*M(1,4)-theta25(29)*M(1,5)-theta25(30)*M(1,2)*M(1,2)-theta25(31)*M(1,2)*M(1,3)-theta25(32)*M(1,3)*M(1,2)-...
    theta25(33)*M(1,3)*M(1,3)-theta25(34)*M(1,2)*M(1,2)-theta25(35)*M(1,2)*M(1,3)-theta25(36)*M(1,3)*M(1,2)-theta25(37)*M(1,3)*M(1,3)-...
    theta25(38)*M(1,2)*M(1,2)-theta25(39)*M(1,2)*M(1,3)-theta25(40)*M(1,3)*M(1,2)-theta25(41)*M(1,3)*M(1,3)-theta25(42)*M(1,4)-theta25(43)*M(1,5)-...
    theta25(44)*M(1,2)*M(1,2)-theta25(45)*M(1,2)*M(1,3)-theta25(46)*M(1,3)*M(1,2)-theta25(47)*M(1,3)*M(1,3)-theta25(48)*M(1,2)*M(1,2)-...
    theta25(49)*M(1,2)*M(1,3)-theta25(50)*M(1,3)*M(1,2)-theta25(51)*M(1,3)*M(1,3)-theta25(52)*M(1,2)*M(1,2)-theta25(53)*M(1,2)*M(1,3)-...
    theta25(54)*M(1,3)*M(1,2)-theta25(55)*M(1,3)*M(1,3)-theta25(56)*M(1,2)*M(1,2)-theta25(57)*M(1,2)*M(1,3)-theta25(58)*M(1,3)*M(1,2)-...
    theta25(59)*M(1,3)*M(1,3)-theta25(60)*M(1,4)-theta25(61)*M(1,5);

% Plug in new coefficients
parfor i = 1:ncoef
    djk  = zeros(ncoef,1);
    dajk = zeros(ncoef,1);
    dbjk = zeros(ncoef,1);
    dcjk = zeros(ncoef,1);
    ddjk = zeros(ncoef,1);
    for k = 1:ncoef
        djk(k)  = theta21(k)*M(index(i,1),indexth(k,1))*M(index(i,2),indexth(k,2))*M(index(i,3),indexth(k,3))*M(index(i,4),indexth(k,4))*M(index(i,5),indexth(k,5));
        dajk(k) = theta22(k)*M(index(i,1),indexth(k,1))*M(index(i,2),indexth(k,2))*M(index(i,3),indexth(k,3))*M(index(i,4),indexth(k,4))*M(index(i,5),indexth(k,5));
        dbjk(k) = theta23(k)*M(index(i,1),indexth(k,1))*M(index(i,2),indexth(k,2))*M(index(i,3),indexth(k,3))*M(index(i,4),indexth(k,4))*M(index(i,5),indexth(k,5));
        dcjk(k) = theta24(k)*M(index(i,1),indexth(k,1))*M(index(i,2),indexth(k,2))*M(index(i,3),indexth(k,3))*M(index(i,4),indexth(k,4))*M(index(i,5),indexth(k,5));
        ddjk(k) = theta25(k)*M(index(i,1),indexth(k,1))*M(index(i,2),indexth(k,2))*M(index(i,3),indexth(k,3))*M(index(i,4),indexth(k,4))*M(index(i,5),indexth(k,5));
    end
    d  = sum(djk);
    da = sum(dajk);
    db = sum(dbjk);
    dc = sum(dcjk);
    dd = sum(ddjk);
    %c,h,b
    ct   = exp(d);
    ht1  = exp(da);
    bt1  = max(db,0)^2;
    u_b  = max(-db,0)^2;
    dht1 = min(max(dc,0)^2,1);
    u_d  = max(-dc,0)^2;
    dat1 = min(max(dd,0)^2,1);
    u_a  = max(-dd,0)^2;
    % Financial wealth
    Vpht1 = rhoph*pht(index(i,3))+(1-rhoph)*log(meanph)+wgridph;
    Vpet1 = rhope*pet(index(i,4))+(1-rhope)*log(meanpe)+wgridpe;
    Vpat1 = rhopa*pat(index(i,5))+(1-rhopa)*log(meanpa)+wgridpa;
    fwt1  = zeros(nstate,nstate,nstate);
    for i_fw = 1:nstate
        for j_fw = 1:nstate
            for k_fw = 1:nstate
                % Next period financial wealth
                fwt1(i_fw,j_fw,k_fw) = log(meanRB*bt1+exp(Vpht1(i_fw))*ht1+ht1*exp(Vpet1(j_fw))*(1-dht1-dat1)+dat1*exp(Vpat1(k_fw))*ht1); 
            end
        end
    end
    fwt1_vec  = fwt1(:);
    root_fwt1 = transfo(fwt1_vec,min(fwt1_vec),max(fwt1_vec)); %vector of [-1,1]
    MFW1      = cheb(root_fwt1,0:npoly); %row = num_rootr column = num_poly
    MFW1      = reshape(MFW1',nnode,nstate,nstate,nstate); %degree*k_fw*j_fw*i_fw
    %log(y_t+1)
    Vyt1      = rhoy*yt(index(i,2))+wgridy; %next period log(y)
    root_yt1  = transfo(Vyt1,min(Vyt1),max(Vyt1)); %vector within [-1,1]
    MY1       = cheb(root_yt1,0:npoly); %row = #rooty column = #poly
    %log(p_t+1^H)
    root_pht1 = transfo(Vpht1,min(Vpht1),max(Vpht1));
    MPH1      = cheb(root_pht1,0:npoly);
    %log(p_t+1^E)
    root_pet1 = transfo(Vpet1,min(Vpet1),max(Vpet1));
    MPE1      = cheb(root_pet1,0:npoly);
    %log(p_t+1^A)
    root_pat1 = transfo(Vpat1,min(Vpat1),max(Vpat1));
    MPA1      = cheb(root_pat1,0:npoly);
    
    %next period control
    VD  = [];
    VDA = [];
    VDC = [];
    djk1  = zeros(ncoef,1);
    dajk1 = zeros(ncoef,1);
    dcjk1 = zeros(ncoef,1);
    for i_state = 1:nstate %ph
        for j_state = 1:nstate %y
            for k_state = 1:nstate %pe
                for l_state = 1:nstate %pa
                    for k = 1:ncoef
                        djk1(k) = theta21(k)*MFW1(indexth(k,1),l_state,k_state,i_state)*MY1(j_state,indexth(k,2))*MPH1(i_state,indexth(k,3))*MPE1(k_state,indexth(k,4))*MPA1(l_state,indexth(k,5));
                        dajk1(k) = theta22(k)*MFW1(indexth(k,1),l_state,k_state,i_state)*MY1(j_state,indexth(k,2))*MPH1(i_state,indexth(k,3))*MPE1(k_state,indexth(k,4))*MPA1(l_state,indexth(k,5));
                        dcjk1(k) = theta24(k)*MFW1(indexth(k,1),l_state,k_state,i_state)*MY1(j_state,indexth(k,2))*MPH1(i_state,indexth(k,3))*MPE1(k_state,indexth(k,4))*MPA1(l_state,indexth(k,5));
                    end
                    d1   = sum(djk1);
                    VD   = [VD;d1];
                    da1  = sum(dajk1);
                    VDA  = [VDA;da1];
                    dc1  = sum(dcjk);
                    VDC  = [VDC;dc1];
                end
            end 
        end
    end
    ct1	 =  exp(VD); % With dimension nstate^2 X 1
    ht2 = exp(VDA);
    dt2 = min(max(VDC,0).^2,1);
    
    % Expected terms
    Et1         = (((ct1.^alpha).*((dt2.*ht2).^(1-alpha))).^(-sigma)).*((dt2.*ht2).^(1-alpha))*alpha.*(ct1.^(alpha-1));
    Et1_reshape = reshape(Et1,nstate,nstate,nstate,nstate); %pe x y x ph
    for i_nstate = 1:nstate
        for j_nstate = 1:nstate
            for k_nstate = 1:nstate
                Et1_reshape(k_nstate,j_nstate,:,i_nstate) = Et1_reshape(k_nstate,j_nstate,:,i_nstate)*(exp(Vpht1(i_nstate))+(1-dht1-dat1)*exp(Vpet1(j_nstate))+dat1*exp(Vpat1(k_nstate)));
            end
        end
    end
    Et1_reshape  = Et1_reshape(:);
    Et1_reshape1 = reshape(Et1,nstate,nstate,nstate,nstate); %pa pe y ph
    for i_nstate = 1:nstate
        Et1_reshape1(:,i_nstate,:,:) = Et1_reshape1(:,i_nstate,:,:)*(exp(Vpet1(i_nstate))*ht1);
    end
    Et1_reshape1 = Et1_reshape1(:);
    
    Et1_reshape2 = reshape(Et1,nstate,nstate,nstate,nstate); %pa pe y ph
    for i_nstate = 1:nstate
        for j_nstate = 1:nstate
            Et1_reshape2(j_nstate,i_nstate,:,:) = Et1_reshape2(j_nstate,i_nstate,:,:)*((exp(Vpat1(j_nstate))-exp(Vpet1(i_nstate)))*ht1);
        end
    end
    Et1_reshape2 = Et1_reshape2(:);
    
    % Euler equation residuals
    rhs1   = beta*meanRB*pi^(-ne/2)*weight'*Et1; %bond
    rhs2   = beta*pi^(-ne/2)*weight'*Et1_reshape; %house
    rhs3   = beta*pi^(-ne/2)*weight'*Et1_reshape1; %deltah
    rhs4   = beta*pi^(-ne/2)*weight'*Et1_reshape2; %deltaa
    rhs5   = exp(fwt(index(i,1)))+exp(yt(index(i,2))); %budget
    lambda =  (ct^alpha*(dht1*ht1)^(1-alpha))^(-sigma)*(dht1*ht1)^(1-alpha)*alpha*ct^(alpha-1);
    lhs2   = lambda*exp(pht(index(i,3)))-(ct^alpha*(dht1*ht1)^(1-alpha))^(-sigma)*ct^alpha*(1-alpha)*(dht1*ht1)^(-alpha)*dht1;
    lhs3   = (ct^alpha*(dht1*ht1)^(1-alpha))^(-sigma)*(dht1*ht1)^(-alpha)*(1-alpha)*ct^(alpha)*ht1-u_d;
    lhs4   = u_a;
    lhs5   =  bt1+ct+ht1*exp(pht(index(i,3)));
    
    LHS1(i)	 = lambda-u_b;
    LHS2(i)	 = lhs2;
    LHS3(i)	 = lhs3;
    LHS4(i)	 = lhs4;
    LHS5(i)	 = lhs5;
    RHS1(i)	 = rhs1;
    RHS2(i)  = rhs2;
    RHS3(i)  = rhs3;
    RHS4(i)  = rhs4;
    RHS5(i)  = rhs5;
end

res1 = LHS1-RHS1;
res2 = LHS2-RHS2;
res3 = LHS3-RHS3;
res4 = LHS4-RHS4;
res5 = LHS5-RHS5;
res  = [res1,res2,res3,res4,res5];
res	 = res(:);
end