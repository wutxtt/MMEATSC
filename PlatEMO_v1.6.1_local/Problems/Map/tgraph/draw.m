global mu req fname
fname = 'dddddd';
global oev1 oevt1 oev2
% [0,0,0]黑色
PopDec(1,1) = 10.73;
PopDec(1,2) = 10.26;
% [0,1,1]青色
% PopDec(1,1) = 2.787;
% PopDec(1,2) = 5.066;
% [0,0,1]蓝色
% PopDec(1,1) = 0.03773;
% PopDec(1,2) = 0.8387;
% PopDec(1,1) = 0.0088;
% PopDec(1,2) = 0.6903;
%[0 0.45 0.74]浅蓝
PopDec(1,1) = 0.004603781942622;
PopDec(1,2) = 6.180859485264221;

dtr = pi / 180.0;
rtd = 180.0 / pi;
mu=398600.4415;
req=6378.1363;
% initial orbit orbital
oev1 = zeros(6,1);
oev1(1) = 6653.14;
oev1(2) = 0.0;
oev1(3) = dtr * 51.6;
oev1(4) = dtr * 0;
oev1(5) = dtr * 0;
oev1(6) = dtr * 96.036 * PopDec(1,1)+pi;
% final orbit orbital
oev2 = zeros(6,1);
oev2(1) = 26553.071184;
oev2(2) = 0.667;
oev2(3) = dtr * 63.4;
oev2(4) = dtr * 270.0;
oev2(5) = dtr * 100;
oev2(6) = dtr * 12.024 * (PopDec(1,1)+PopDec(1,2))+2*pi;
% transfer orbit - first impulse
[ri1, vi1] = orb2eci(mu, oev1);
[ri, vi] = orb2eci(mu, oev1);
[rf, vf] = orb2eci(mu, oev2);
for j = 1:1:3
    sv1(j) = ri(j);
    sv1(j + 3) = vi(j);
    sv2(j) = rf(j);
    sv2(j + 3) = vf(j);
end
dtof = PopDec(1,2)*60.0*60.0;
nrev = 0;
[vito, vfto] = glambert(mu, sv1, sv2, dtof, nrev);
oevt1 = eci2orb1 (mu, ri1, vito');
tgraphics(1,dtof);
