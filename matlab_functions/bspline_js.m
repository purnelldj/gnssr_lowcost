function [res] = bspline_js(coefs_0,t_all,sinelv_all,snr_all,knots,bspline_order,...
    satno_all,gps,glo,gal,antno_all,meanhgts,sig)

t_alls=t_all;
sinelv_alls=sinelv_all;
snr_alls=snr_all;
satno_alls=satno_all;

res=[];
ants=unique(antno_all);

consts=gps+glo+gal;
height_coeff=coefs_0(1:end-consts*2-1);

for k=1:numel(ants)

now=antno_all(:)==k;
t_all=t_alls(now);
sinelv_all=sinelv_alls(now);
snr_all=snr_alls(now);
satno_all=satno_alls(now);

h1=bspline_deboor(bspline_order+1,knots,height_coeff,t_all);
h1=h1.';
if k>1
    h1=h1+meanhgts(k)-meanhgts(1);
end

tmpc=0;
% GPS
if gps==1
tmpc=tmpc+1;
gps=satno_all(:,1)<33;
if sig==1
Lcar=(299792458/(1575.42e06)); % for GPS
elseif sig==2
Lcar=299792458/1227.60e06;
end
Lk=(2*pi)/Lcar;
modelSNR = (coefs_0(end-1)*sin((4*pi*h1(gps).*sinelv_all(gps))/Lcar)+...
    coefs_0(end-2)*cos((4*pi*h1(gps).*sinelv_all(gps))/Lcar)).*...
    exp(-4*Lk^2*coefs_0(end)*sinelv_all(gps).^2);
res=[res;modelSNR-snr_all(gps)];
end
% GLO
if glo==1
tmpc=tmpc+1;
glo=satno_all(:,1)>32 & satno_all(:,1)<57;
load('glonasswlen.mat')
satnotmp=satno_all(glo,1);
for ij=1:numel(satnotmp)
    if sig==1
    Lcar(ij)=glonassl1(satnotmp(ij)-32);
    elseif sig==2
    Lcar(ij)=glonassl2(satnotmp(ij)-32);
    end
end
Lcar=Lcar.';
Lk=(2*pi)./Lcar;
modelSNR = (coefs_0(end-(tmpc-1)*2-1)*sin((4*pi*h1(glo).*sinelv_all(glo))./Lcar)+...
    coefs_0(end-(tmpc-1)*2-2)*cos((4*pi*h1(glo).*sinelv_all(glo))./Lcar)).*...
    exp(-4*Lk.^2*coefs_0(end).*sinelv_all(glo).^2);
res=[res;modelSNR-snr_all(glo)];
end
% GAL
if gal==1
tmpc=tmpc+1;
gal=satno_all(:,1)>56;
if sig==1
Lcar=299792458/1575.42e06;
elseif sig==2
Lcar=299792458/1227.60e06;
end
Lk=(2*pi)/Lcar;
modelSNR = (coefs_0(end-(tmpc-1)*2-1)*sin((4*pi*h1(gal).*sinelv_all(gal))/Lcar)+...
    coefs_0(end-(tmpc-1)*2-2)*cos((4*pi*h1(gal).*sinelv_all(gal))/Lcar)).*...
    exp(-4*Lk^2*coefs_0(end)*sinelv_all(gal).^2);
res=[res;modelSNR-snr_all(gal)];
end
end

end

