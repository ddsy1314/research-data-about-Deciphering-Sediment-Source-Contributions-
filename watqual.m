thgra = 1.047;thrho = 1.047;thrs1 = 1.024;thrs2 = 1.074;thrs3 = 1.074;
thrs4 = 1.024;thrs5 = 1.024;thbc1 = 1.083;thbc2 = 1.047;thbc3 = 1.047;thbc4 = 1.047;
thrk1 = 1.047;thrk2 = 1.024;thrk3 = 1.024;thrk4 = 1.060;
jrch=inum1;
dcoef=3;
wtrin=inputvalue(3);
if wtrin>10^-4
    chlin = 0;
    algin = 0;
    orgnin = 0;
    ammoin = 0;
    nitritin = 0;
    nitratin = 0;
    orgpin = 0;
    dispin = 0;
    cbodin = 0;
    disoxin = 0;
    n1=0;n2=0;n3=0;n4=0;n5=0;n6=0;n7=0;n8=0;
    p1=0;p2=0;p3=0;p4=0;p5=0;
    cinn = 0;
    if wtrin > 0.001
        chlin = 1000*inputvalue(6)/ wtrin;
        algin = 1000* chlin / ai0;
        orgnin = 1000* inputvalue(9)/ wtrin;
        ammoin = 1000* inputvalue(12)/ wtrin;
        nitritin = 1000* inputvalue(5)/ wtrin;
        nitratin = 1000* inputvalue(11)/ wtrin;
        orgpin = 1000* inputvalue(10)/ wtrin;
        dispin = 1000* inputvalue(13)/ wtrin;
        cbodin = 1000* inputvalue(7)/ wtrin;
        disoxin = 1000*inputvalue(8)/ wtrin;
    end
    wtrtot = 0;
    algcon = 0;
    orgncon = 0;
    nh3con = 0;
    no2con = 0;
    no3con = 0;
    orgpcon = 0;
    solpcon = 0;
    cbodcon = 0;
    o2con = 0;
    rch_cbod=max(10^-6,rch_cbod);
    wtrtot=wtrin+rchwtr;
    algcon = (algin * wtrin + algae(inum1) * rchwtr) / wtrtot;
    orgncon = (orgnin * wtrin + organicn (inum1)* rchwtr) / wtrtot;
    nh3con = (ammoin * wtrin + ammonian(inum1) * rchwtr) / wtrtot;
    no2con = (nitritin * wtrin + nitriten(inum1) * rchwtr) / wtrtot;
    no3con = (nitratin * wtrin + nitraten(inum1) * rchwtr) / wtrtot;
    orgpcon = (orgpin * wtrin + organicp(inum1) * rchwtr) / wtrtot;
    solpcon = (dispin * wtrin + disolvp(inum1) * rchwtr) / wtrtot;
    cbodcon = (cbodin * wtrin + rch_cbod(inum1) * rchwtr) / wtrtot;
    o2con = (disoxin * wtrin + rch_dox(inum1) * rchwtr) / wtrtot;
    if orgncon<10^-6
        orgncon=0;
    end
    if nh3con<10^-6
        nh3con=0;
    end
    if no2con<10^-6
        no2con=0;
    end
    if no3con<10^-6
        no3con=0;
    end
    if orgpcon<10^-6
        orgpcon=0;
    end
    if solpcon<10^-6
        solpcon=0;
    end
    if cbodcon<10^-6
        cbodcon=0;
    end
    if o2con<10^-6
        o2con=0;
    end
    wtmp=0;
    wtmp=5+0.75*inputvalue(14);
    if wtmp<=0
        wtmp=0.1;
    end
    cinn=nh3con+no3con;
    ww=0;
    xx=0;
    yy=0;
    zz=0;
    ww=-139.34410+(1.575701*10^5/(wtmp+273.15));
    xx=6.642308*10^7/((wtmp+273.15)^2);
    yy=1.243800*10^10/((wtmp+273.15)^3);
    zz=8.621949*10^11/((wtmp+273.15)^4);
    soxy=exp(ww-xx+yy-zz);
    if soxy<10^-6
        soxy=0;
    end
    cordo=0;
    if o2con<=0.001
        o2con<=0.001;
    end
    if o2con>30
        o2con=30;
    end
    cordo=1-exp(-0.6*o2con);
    bc1mod=0;
    bc2mod=0;
    bc1mod=bc1*cordo;
    bc2mod=bc2*cordo;
    tday=0;
    tday=rttime/24;
    if tday>1
        tday=1;
    end
    if ai0*algcon>10^-6
        lambda=lambda0+(lambda1*ai0*algcon)+lambda2*(ai0*algcon)^0.66667;
    else
        lambda=lambda0;
    end
    if lambda>lambda0
        lambda=lambda0;
    end
    fnn=0;
    fpp=0;
    fnn=cinn/(cinn+k_n);
    fpp=solpcon/(solpcon+k_p);
    algi=0;
    run('clgen.m')
    if dayl>0
        algi=inputvalue(16)*tfact/dayl;
    else
        algi=0;
    end
    fl_1=0;
    fl1=0;
    fl_1=(1/(lambda*rchdep))*log((k_l+algi)/(k_l+algi*(exp(-lambda*rchdep))));
    fl1=0.92*(dayl/24)*fl_1;
    gra=0;
    if igropt==1
        gra=mumax*fl1*fnn*fpp;
    elseif igropt==2
        gra=mumax*fl1*min(fnn,fpp);
    elseif igropt==3
        if fnn>10^-6 && fpp>10^-6
            gra=mumax*fl1*2/((1/fnn)+(1/fpp));
        else
            gra=0;
        end
    end
    algae(inum1)=0;
    algae(inum1)=algcon+(theta(gra,thgra,wtmp)*algcon-theta(rhoq,thrho,wtmp)*algcon-theta(rs1,thrs1,wtmp)/rchdep*algcon)*tday;
    if algae(inum1)<10^-6
        algae(inum1)=0;
    end
    if algae(inum1)>5000
        algae(inum1)=5000;
    end
    if algae(inum1)>dcoef*algcon
        algae(inum1)=dcoef*algcon;
    end
    chlora(inum1)=0;
    chlora(inum1)=algae(inum1)*ai0/1000;
    yy=0;
    zz=0;
    yy=theta(rk1,thrk1,wtmp)*cbodcon;
    zz=theta(rk3,thrk3,wtmp)*cbodcon;
    rch_cbod(inum1)=0;
    rch_cbod(inum1)=cbodcon-(yy+zz)*tday;
    coef=exp(-theta(rk1,thrk1,wtmp)*tday);
    cbodrch=coef*cbodcon;
    coef=exp(-theta(rk3,thrk3,wtmp)*tday);
    cbodrch=coef*cbodrch;
    rch_cbod(inum1)=cbodrch;
    if rch_cbod(inum1)<10^-6
        rch_cbod(inum1)=0;
    end
    if rch_cbod(inum1)>dcoef*cbodcon
        rch_cbod(inum1)=dcoef*cbodcon;
    end
    uu=0;
    vv = 0;
    ww = 0;
    xx = 0;
    yy = 0;
    zz = 0;
    rhoq = 1.0;
    rk2 = 1.0;
    uu=theta(rk2,thrk2,wtmp)*(soxy-o2con);
    vv=(ai3*theta(gra,thgra,wtmp)-ai4*theta(rhoq,thrho,wtmp))*algcon;
    ww=theta(rk1,thrk1,wtmp)*cbodcon;
    xx=theta(rk4,thrk4,wtmp)/(rchdep*1000);
    yy=ai5*theta(bc1mod,thbc1,wtmp)*nh3con;
    zz=ai6*theta(bc2mod,thbc2,wtmp)*no2con;
    rch_dox(inum1)=o2con+(uu+vv-ww-xx-yy-zz)*tday;
    rch_dox(inum1)=min(0.1,rch_dox(inum1));
    doxrch=soxy;
    coef=exp(-0.1*ww);
    doxrch=coef*doxrch;
    coef=1-(theta(rk4,thrk4,wtmp)/100);
    doxrch=coef*doxrch;
    coef=exp(-0.05*yy);
    doxrch=coef*doxrch;
    coef=exp(-0.05*zz);
    doxrch=coef*doxrch;
    uu=theta(rk2,thrk2,wtmp)/100*(soxy-doxrch);
    rch_dox(inum1)=doxrch+uu;
    if rch_dox(inum1)<10^-6
        rch_dox(inum1)=0;
    end
    if rch_dox(inum1)>soxy
       rch_dox(inum1)=soxy;
    end
    xx=0;
    yy=0;
    zz=0;
    xx=ai1*theta(rhoq,thrho,wtmp)*algcon;
    yy=theta(bc3,thbc3,wtmp)*orgncon;
    zz=theta(rs4,thrs4,wtmp)*orgncon;
    n7=xx;
    n8=zz;
    organicn(inum1)=0;
    organicn(inum1)=orgncon+(xx-yy-zz)*tday;
    if organicn(inum1)<10^-6
        organicn(inum1)=0;
    end
    if organicn(inum1)>dcoef*orgncon
       organicn(inum1)=dcoef*orgncon;
    end
    f1=0;
    f1=p_n*nh3con/(p_n*nh3con+(1-p_n)*no3con+10^-6);
     ww = 0;
    xx = 0;
    yy = 0;
    zz = 0;
    ww=theta(bc3,thbc3,wtmp)*orgncon;
    xx=theta(bc1mod,thbc1,wtmp)*nh3con;
    yy=theta(rs3,thrs3,wtmp)/(rchdep*1000);
    zz=f1*ai1*algcon*theta(gra,thgra,wtmp);
    n1=ww;
    n2=zz;
    n3=yy;
    n4=xx;
    ammonian(inum1)=0;
    ammonian(inum1)=nh3con+(ww-xx+yy-zz)*tday;
    if ammonian(inum1)<10^-6
        ammonian(inum1)=0;
    end
    if ammonian(inum1)>dcoef*nh3con && nh3con>0
       ammonian(inum1)=dcoef*nh3con;
    end
    yy=0;
    zz=0;
    yy=theta(bc1mod,thbc1,wtmp)*nh3con;
    zz=theta(bc2mod,thbc2,wtmp)*no2con;
    nitriten(inum1)=0;
    nitriten(inum1)=no2con+(yy-zz)*tday;
    if nitriten(inum1)<10^-6
        nitriten(inum1)=0;
    end
    if nitriten(inum1)>dcoef*no2con && no2con>0
       nitriten(inum1)=dcoef*no2con;
    end
    yy=0;
    zz=0;
    yy=theta(bc2mod,thbc2,wtmp)*no2con;
    zz=(1-f1)*ai1*algcon*theta(gra,thgra,wtmp);
    nitraten(inum1)=0;
    n5=yy;
    n6=zz;
    nitraten(inum1)=no3con+(yy-zz)*tday;
    if nitraten(inum1)<10^-6
        nitraten(inum1)=0;
    end
    if nitraten(inum1)>dcoef*no3con && no3con>0
       nitraten(inum1)=dcoef*no3con;
    end
    xx=0;
    yy=0;
    zz=0;
    xx=ai2*theta(rhoq,thrho,wtmp)*algcon;
    yy=theta(bc4,thbc4,wtmp)*orgpcon;
    zz=theta(rs5,thrs5,wtmp)*orgpcon;
    organicp(inum1)=0;
    p1=xx;
    p2=zz;
    p3=yy;
    organicp(inum1)=orgpcon+(xx-yy-zz)*tday;
    if organicp(inum1)<10^-6
        organicp(inum1)=0;
    end
    if organicp(inum1)>dcoef*orgpcon
       organicp(inum1)=dcoef*orgpcon;
    end
    xx=0;
    yy=0;
    zz=0;
    xx=theta(bc4,thbc4,wtmp)*orgpcon;
    yy=theta(rs2,thrs2,wtmp)/(rchdep*1000);
    zz=ai2*theta(rhoq,thrho,wtmp)*algcon;
    p4=yy;
    p5=zz;
    disolvp(inum1)=0;
    disolvp(inum1)=solpcon+(xx+yy-zz)*tday;
    if disolvp(inum1)<10^-6
        disolvp(inum1)=0;
    end
    if disolvp(inum1)>dcoef*solpcon
       disolvp(inum1)=dcoef*solpcon;
    end
else
    algin = 0.0;
        chlin = 0.0;
        orgnin = 0.0;
        ammoin = 0.0;
        nitritin = 0.0;
        nitratin = 0.0;
        orgpin = 0.0;
        dispin = 0.0;
        cbodin = 0.0;
        disoxin = 0.0;
        algae(inum1) = 0.0;
        chlora(inum1) = 0.0;
        organicn (inum1)= 0.0;
        ammonian (inum1)= 0.0;
        nitriten (inum1)= 0.0;
        nitraten (inum1)= 0.0;
        organicp(inum1)= 0.0;
        disolvp(inum1)= 0.0;
        rch_cbod(inum1)= 0.0;
        rch_dox(inum1)= 0.0;
        soxy = 0.0;
        orgncon = 0.0;
        no2con=0;
end
%%
function kk=theta(a,b,c)
kk=0;
kk=a*b^(c-20);
end