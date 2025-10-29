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