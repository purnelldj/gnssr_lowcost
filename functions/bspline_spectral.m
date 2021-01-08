function bsplinenadj = bspline_spectral(coefs,bspline_order,knots,tanthter,t_init,dohgtcor)

% made by david purnell (2020)

evendt=1/(24*60);
t_even=min(knots):evendt:max(knots); % 5 min intervals, can change

bsplineeven=bspline_deboor(bspline_order+1,knots,coefs,t_even);

dhdteven=gradient(bsplineeven,evendt); % currently in m/day so
dhdt=interp1(t_even,dhdteven,t_init,'linear'); % can dry different techniques
dhdt=dhdt./86400; % now in m/s (?)
dhdt=dhdt.';
tanthter=tanthter.';

if dohgtcor==1
bsplinenadj=bspline_deboor(bspline_order+1,knots,coefs,t_init)+dhdt.*tanthter;
else
bsplinenadj=bspline_deboor(bspline_order+1,knots,coefs,t_init);
end

end

