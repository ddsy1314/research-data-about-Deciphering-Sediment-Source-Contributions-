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