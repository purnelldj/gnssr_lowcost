function bsplinenadj = bspline_spectral_accel(coefs,bspline_order,knots,tanthter,xinit,tinit,dohgtcor)

% written by Dave Purnell https://github.com/purnelldj/gnssr_lowcost

evendt=1/(24*60);
t_even=min(knots):evendt:max(knots); % 5 min intervals, can change

bsplineeven=bspline_deboor(bspline_order+1,knots,coefs,t_even);

dhdteven=gradient(bsplineeven,evendt); % currently in m/day so
dhdteven=dhdteven./86400; % now in m/s
dhdt=interp1(t_even,dhdteven,xinit,'linear'); % can dry different techniques
dhdt=dhdt.';
dhdt2even=gradient(dhdteven,evendt*86400);
dhdt2=interp1(t_even,dhdt2even,xinit,'linear'); % should be m/s^2
dhdt2=dhdt2.';
tanthter=tanthter.';
tinit=tinit*86400; % turn into seconds

if dohgtcor==1
bsplinenadj=bspline_deboor(bspline_order+1,knots,coefs,xinit)+(dhdt+tinit.*dhdt2).*tanthter;
else
bsplinenadj=bspline_deboor(bspline_order+1,knots,coefs,xinit);
end

end

