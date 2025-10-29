result=[];
sub=1
% percent=1
parfor sub=1:48
     for percent=1:6
        per=1-(percent-1)*0.1;
        rchstor=zeros(1,48);
        bankst=zeros(1,48);
        sedst=zeros(1,48);
        ch_orgn=zeros(1,48);
        ch_orgp=zeros(1,48);
        rch_cbod=zeros(1,48);
        algae=zeros(1,48);
        organicn=zeros(1,48);
        nitriten=zeros(1,48);
        ammonian=zeros(1,48);
        nitraten=zeros(1,48);
        organicp=zeros(1,48);
        disolvp=zeros(1,48);
        rch_dox=zeros(1,48);
        depch=zeros(1,48);
        chlora=zeros(1,48);
        out=[];
        %%
        
        day=1;
        date=datetime(1996,12,31)
        while date<=datetime(2020,12,30)
            % while day<=31
            date=datetime(1996,12,31)+day;
            input_values=zeros(143,16);
            input_values(1:48,:)=inputs((day-1)*48+1:(day-1)*48+48,:);
            for i=1:1:size(structure,1)
                inum1=structure(i,3);
                inum2=structure(i,4);
                if structure(i,1)==2
                    inputvalue=[inum1 inputs((day-1)*48+inum1,2) input_values(inum2,3:13) inputs((day-1)*48+inum1,14:16)];
                    if inum1==sub
                        inputvalue(3)=per*inputvalue(3);
                    end
                    chd=paratemer_rte(inum1,1);chn2=paratemer_rte(inum1,2);chs2=paratemer_rte(inum1,3);chw2=paratemer_rte(inum1,4);chl2=paratemer_rte(inum1,5);
                    chk2=paratemer_rte(inum1,6);evrch=paramater_bsn(1);trnsrch=paramater_bsn(2);sub_ha=paratemer_rte(inum1,8)*100;ch_revap=paratemer_rte(inum1,9);
                    alaph_bnk=paratemer_rte(inum1,10);alaph_bnke=exp(-alaph_bnk);prf=paramater_bsn(21);spcon=paramater_bsn(22);spexp=paramater_bsn(23);
                    ch_cov2=paratemer_rte(inum1,7);ch_onco=paratemer_rte(inum1,11);ch_opco=paratemer_rte(inum1,12);bc1=paratemer_rte(inum1,13);
                    bc2=paratemer_rte(inum1,14);ai0=paramater_bsn(3);lambda0 = paramater_bsn(4);lambda1 =paramater_bsn(5);lambda2 = paramater_bsn(6);
                    k_n=paramater_bsn(7);k_p=paramater_bsn(8);latitude=paratemer_rte(inum1,15);tfact=paramater_bsn(9);k_l=paramater_bsn(10);
                    igropt=paramater_bsn(11);mumax=paramater_bsn(12);rhoq=paramater_bsn(13);rs1=paratemer_rte(inum1,16);rk1=paratemer_rte(inum1,17);
                    rk3=paratemer_rte(inum1,18);rk4=paratemer_rte(inum1,19);ai1=paramater_bsn(14);ai2=paramater_bsn(15);ai3=paramater_bsn(16);
                    ai4=paramater_bsn(17);ai5=paramater_bsn(18);ai6=paramater_bsn(19);bc3=paratemer_rte(inum1,20);bc4=paratemer_rte(inum1,21);
                    rs2=paratemer_rte(inum1,22);rs3=paratemer_rte(inum1,23);rs4=paratemer_rte(inum1,24);rs5=paratemer_rte(inum1,25);
                    p_n=paramater_bsn(20);ch_side=paramater_bsn(24);
                    rchwtr = rchstor(inum1);
                    sedrch=0;
                    pteday=inputvalue(2);
                    wtrin=inputvalue(3);
                    vol=wtrin+rchstor(inum1);
                    volrt=vol/84600;
                    c=ch_side;
                    b=chw2-2*chd*c;
                    phi6=b;
                    p=phi6+2*chd*sqrt(1+c^2);
                    phi1=b+chd+c*chd^2;
                    rh=phi1/p;
                    maxrt=Qman(phi1,rh,chn2,chs2);
                    sdti=0;
                    rchdep=0;
                    p=0;
                    rh=0;
                    vc=0;
                    if volrt>maxrt
                        rcharea=phi1;
                        rchdep=chd;
                        p=phi6+2*chd*sqrt(1+c*c);
                        rh=phi1/p;
                        sdti=maxrt;
                        adddep=0;
                        while sdti<volrt
                            adddep=adddep+0.01;
                            addarea=rcharea+(chw2*5+4*adddep)*adddep;
                            addp=p+(chw2*4)+2*adddep*sqrt(1+4*4);
                            rh=addarea/addp;
                            sdti=Qman(addarea,rh,chn2,chs2);
                        end
                        rcharea=addarea;
                        rchdep=chd+adddep;
                        p=addp;
                        sdti=volrt;
                    else
                        while sdti<volrt
                            rchdep=rchdep+0.01;
                            rcharea=(phi6+c*rchdep)*rchdep;
                            p=phi6+2*rchdep*sqrt(1+c*c);
                            rh=rcharea/p;
                            sdti=Qman(rcharea,rh,chn2,chs2);
                        end
                        sdti=volrt;
                    end
                    topw=0;
                    if rchdep<=chd
                        topw=phi6+2*rchdep*c;
                    else
                        topw=5*chw2+2*(rchdep-chd)*4;
                    end
                    det=24;
                    if sdti>0
                        vc=sdti/rcharea;
                        velchan=vc;
                        rttime=chl2*1000/(3600*vc);
                        scoef=0;
                        rtwtr=0;
                        scoef=det/(rttime+det);
                        if scoef>1
                            scoef=1;
                        end
                        rtwtr=scoef*(wtrin+rchstor(inum1));
                        rchstor(inum1)=rchstor(inum1)+wtrin-rtwtr;
                        if rchstor<0
                            rchstor=0;
                        end
                        rttlc=0;
                        if rtwtr>0
                            rttlc=det*chk2*chl2*p;
                            rttlc2=rttlc*rchstor(inum1)/(rtwtr+rchstor(inum1));
                            if rchstor(inum1)<rttlc2
                                rttlc2=min(rttlc2,rchstor(inum1));
                                rchstor(inum1)=rchstor(inum1)-rttlc2;
                                rttlc1=rttlc-rttlc2;
                                if rtwtr<rttlc1
                                    rttlc1=min(rttlc1,rtwtr);
                                    rtwtr=rtwtr-rttlc1;
                                else
                                    rtwtr=rtwtr-rttlc1;
                                end
                            else
                                rchstor(inum1)=rchstor(inum1)-rttlc2;
                                rttlc1=rttlc-rttlc2;
                                if rtwtr<=rttlc1
                                    rttlc1=min(rttlc1,rtwtr);
                                    rtwtr=rtwtr-rttlc1;
                                else
                                    rtwtr=rtwtr-rttlc1;
                                end
                            end
                            rttlc=rttlc1+rttlc2;
                        end
                        rtevp=0;
                        if (rtwtr>0)
                            aaa=evrch*inputvalue(2)/1000;
                            if rchdep<chd
                                rtevp=aaa*chl2*1000*topw;
                            else
                                if aaa<(rchdep-chd)
                                    rtevp=aaa*chl2*1000*topw;
                                else
                                    rtevp=rchdep-chd;
                                    rtevp=rtevp+(aaa-(rchdep-chd));
                                    topw=phi6+2*chd*c;
                                    rtevp=rtevp*chl2*1000*topw;
                                end
                            end
                            rtevp2=rtevp*rchstor/(rtwtr+rchstor);
                            if rchstor(inum1)<rtevp2
                                rtevp2=min(rtevp2,rchstor(inum1));
                                rchstor(inum1)=rchstor(inum1)-rtevp2;
                                rtevp1=rtevp-rtevp2;
                                if rtwtr<=rtevp1
                                    rtevp1=min(rtevp1,rtwtr);
                                    rtwtr=rtwtr-rtevp1;
                                else
                                    rtwtr=rtwtr-rtevp1;
                                end
                            else
                                rchstor(inum1)=rchstor(inum1)-rtevp2;
                                rtevp1=rtevp-rtevp2;
                                if rtwtr<=rtevp1
                                    rtevp1=min(rtevp1,rtwtr);
                                    rtwtr=rtwtr-rtevp1;
                                else
                                    rtwtr=rtwtr-rtevp1;
                                end
                            end
                            rtevp=rtevp1+rtevp2;
                        end
                    else
                        rtwtr=0;
                        sdti=0;
                        rchstor(inum1)=0;
                        velchan=0;
                        flwin=0;
                        flwout=0;
                    end
                    volinprev=wtrin;
                    qoutprev=rtwtr;
                    if rtwtr<0
                        rtwtr=0;
                    end
                    if rchstor(inum1)<0
                        rchstor(inum1)=0;
                    end
                    if rchstor(inum1)<10
                        rtwtr=rtwtr+rchstor(inum1);
                        rchstor(inum1)=0;
                    end
                    
                    if rttlc>0
                        bankst(inum1)=bankst(inum1)+rttlc*(1-trnsrch);
                        if sub_ha>10^-9
                            subwtr=rttlc*trnsrch/(sub_ha*10);
                        end
                    end
                    revapday=ch_revap*inputvalue(2)*chl2*chw2;
                    revapday=min(revapday,bankst(inum1));
                    bankst(inum1)=bankst(inum1)-revapday;
                    qdbank=bankst(inum1)*(1-alaph_bnke);
                    bankst(inum1)=bankst(inum1)-qdbank;
                    rtwtr=rtwtr+qdbank;
                    if inum1==inum2
                        if rtwtr>0 && rchdep>0
                            sedrch=inputvalue(4);
                        end
                    else
                        sedin=0;
                        if rtwtr>0 && rchdep>0
                            qdin=0;
                            qdin=rtwtr+rchstor(inum1);
                            if qdin>0.01
                                sedin=0;
                                sedin=inputvalue(4)+sedst(inum1);
                                sedinorg=sedin;
                                peakr=prf*sdti;
                                vc=0;
                                if rchdep<0.01
                                    vc=0.01;
                                else
                                    vc=peakr/rcharea;
                                end
                                if vc>5
                                    vc=5;
                                end
                                tbase=chl2*1000/(3600*24*vc);
                                if tbase>1
                                    tbase=1;
                                end
                                cyin=0;
                                cych=0;
                                depnet=0;
                                deg=0;
                                deg1=0;
                                deg2=0;
                                dep=0;
                                cyin=sedin/qdin;
                                cych=spcon*vc^spexp;
                                depnet=rtwtr*(cych-cyin);
                                if abs(depnet)<1*10^-6
                                    depnet=0;
                                end
                                if depnet>10^-6
                                    deg=depnet;
                                    if deg>=depch(inum1)
                                        deg1=depch(inum1);
                                        deg2=(deg-deg1)*ch_erodmo(inum1,month(date))*ch_cov2;
                                    else
                                        deg1=deg;
                                        deg2=0;
                                    end
                                    dep=0;
                                else
                                    dep=-depnet;
                                    deg=0;
                                    deg1=0;
                                    deg2=0;
                                end
                                depch(inum1)=depch(inum1)+dep-deg1;
                                if depch(inum1)<10^-6
                                    depch(inum1)=0;
                                end
                                sedin=sedin+deg1+deg2-dep;
                                if sedin<10^-6
                                    sedin=0;
                                end
                                outfract=rtwtr/qdin;
                                if outfract>1
                                    outfract=1;
                                end
                                sedrch=sedin*outfract;
                                if sedrch<10^-6
                                    sedrch=0;
                                end
                                sedst(inum1)=sedin-sedrch;
                                if sedst(inum1)<10^-6
                                    sedst(inum1)=0;
                                end
                                ch_orgn(inum1)=deg2*ch_onco/1000;
                                ch_orgp(inum1)=deg2*ch_opco/1000;
                            else
                                sedrch=0;
                                sedst(inum1)=sedin;
                            end
                        end
                    end
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
                        daii=-(datenum(datetime(year(date),1,1))-datenum(date))+1;
                        if inputvalue(15)>=0.1
                            npcp=2;
                        else
                            npcp=1;
                        end
                        sd=0;
                        sd=asin(0.4*sin((daii-82)/58.09));
                        dd=0;
                        dd=1+0.033*cos(daii/58.09);
                        ch=0;
                        h=0;
                        ch=-sind(latitude)*tan(sd)/cosd(latitude);
                        if ch>1
                            h=0;
                        elseif ch>=-1
                            h=acos(ch);
                        else
                            h=3.1416;
                        end
                        dayl=7.6394*h;
                        ys=0;
                        yc=0;
                        ys=sind(latitude)*sin(sd);
                        yc=cosd(latitude)*cos(sd);
                        hru_rmx=30*dd*(h*ys+yc*sin(h));
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
                    orgn=organicn(inum1) *rtwtr/1000+ch_orgn(inum1) ;
                    orgp=organicp(inum1) *rtwtr/1000+ch_orgp(inum1) ;
                    solp=disolvp(inum1) *rtwtr/1000;
                    nh4=ammonian(inum1) *rtwtr/1000;
                    no2=nitriten(inum1) *rtwtr/1000;
                    no3=nitraten(inum1) *rtwtr/1000;
                    do=rch_dox(inum1) *rtwtr/1000;
                    bod=rch_cbod(inum1) *rtwtr/1000;
                    ch_a=chlora(inum1) *rtwtr/1000;
                    out_value=[inum1 inputvalue(2) rtwtr sedrch no2 ch_a bod do orgn orgp no3 nh4 solp inputvalue(14:16)];
                    input_values(structure(i,2),:)=out_value;
                    if inum1==48
                        out=[out;day out_value];
                    end
                elseif structure(i,1)==5
                    input_values(structure(i,2),:)=[input_values(inum1,:)+input_values(inum2,:)];
                end
            end
            day=day+1;
        end
        result(sub,percent)= mean(out( :,5));
     end
end

function k=Qman(a,b,c,d)
k=a*b^0.6666*sqrt(d)/c;
end

function kk=theta(a,b,c)
kk=0;
kk=a*b^(c-20);
end
