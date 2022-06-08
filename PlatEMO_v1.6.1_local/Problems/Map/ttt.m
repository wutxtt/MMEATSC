oev1 = [15000,0,0,0,0,pi];
oev2 = [15000,0,0,0,0,2*pi];
mu=42828;
[ri, vi] = orb2eci(mu, oev1);
[rf, vf] = orb2eci(mu, oev2);
for i = 1:1:3
    
    sv1(i) = ri(i);
    
    sv1(i + 3) = vi(i);

    sv2(i) = rf(i);
    
    sv2(i + 3) = vf(i);
    
end
tof = 3600;
nrev = 0;
[vito, vfto] = glambert(mu, sv1, sv2, tof, nrev);

dvi(1) = vito(1) - vi(1);
dvi(2) = vito(2) - vi(2);
dvi(3) = vito(3) - vi(3);

dvf(1) = vf(1) - vfto(1);
dvf(2) = vf(2) - vfto(2);
dvf(3) = vf(3) - vfto(3);

delta_v = norm(dvi) + norm(dvf);
delta_t = tof;