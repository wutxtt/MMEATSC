function tgraphics(soln, dtof)

% create graphics display and disk file of trajectory graphics

% input

%  soln = solution number
%  dtof = time-of-flight (seconds)

% Orbital Mechanics with MATLAB

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global mu req fname

global oev1 oevt1 oev2

rp1_x = zeros(301, 1);

rp1_y = zeros(301, 1);

rp1_z = zeros(301, 1);

rp2_x = zeros(301, 1);

rp2_y = zeros(301, 1);

rp2_z = zeros(301, 1);

rp3_x = zeros(301, 1);

rp3_y = zeros(301, 1);

rp3_z = zeros(301, 1);

% clear graphics

% clf;

% create disk filename

gfilename1 = horzcat(fname, '_traj', num2str(soln), '.tif');

gfilename2 = horzcat(fname, '_traj', num2str(soln), '.fig');

% compute initial state vectors

[ri1, vi1] = orb2eci(mu, oev1);

[ri2, vi2] = orb2eci(mu, oevt1);

[ri3, vi3] = orb2eci(mu, oev2);

% compute orbital periods

period1 = 2.0 * pi * oev1(1) * sqrt(oev1(1) / mu);

period3 = 2.0 * pi * oev2(1) * sqrt(oev2(1) / mu);

deltat1 = period1 / 300;

simtime1 = -deltat1;

deltat2 = dtof / 300;

simtime2 = -deltat2;

deltat3 = period3 / 300;

simtime3 = -deltat3;

% compute graphics data

for i = 1:1:301
    
    simtime1 = simtime1 + deltat1;
    
    simtime2 = simtime2 + deltat2;
    
    simtime3 = simtime3 + deltat3;
    
    % park orbit "normalized" position vector
    
    [rwrk, ~] = twobody2 (mu, simtime1, ri1, vi1);
    
    rp1_x(i) = rwrk(1) / req;
    
    rp1_y(i) = rwrk(2) / req;
    
    rp1_z(i) = rwrk(3) / req;
    
    % transfer orbit position vector
    
    [rwrk, ~] = twobody2 (mu, simtime2, ri2, vi2);
    
    rp2_x(i) = rwrk(1) / req;
    
    rp2_y(i) = rwrk(2) / req;
    
    rp2_z(i) = rwrk(3) / req;
    
    % final orbit position vector
    
    [rwrk, ~] = twobody2 (mu, simtime3, ri3, vi3);
    
    rp3_x(i) = rwrk(1) / req;
    
    rp3_y(i) = rwrk(2) / req;
    
    rp3_z(i) = rwrk(3) / req;
    
end

% create axes vectors

xaxisx = [1 1.5];
xaxisy = [0 0];
xaxisz = [0 0];

yaxisx = [0 0];
yaxisy = [1 1.5];
yaxisz = [0 0];

zaxisx = [0 0];
zaxisy = [0 0];
zaxisz = [1 1.5];

hold on;

% plot earth

[x, y, z] = sphere(24);

h = surf(x, y, z);

colormap([127/255 1 222/255]);

set (h, 'edgecolor', [1 1 1]);

% plot coordinate system axes

plot3(xaxisx, xaxisy, xaxisz, '-g', 'LineWidth', 1);

plot3(yaxisx, yaxisy, yaxisz, '-r', 'LineWidth', 1);

plot3(zaxisx, zaxisy, zaxisz, '-b', 'LineWidth', 1);

% plot initial orbit

plot3(rp1_x, rp1_y, rp1_z, '-r', 'LineWidth', 1.5);

% plot transfer orbit

plot3(rp2_x, rp2_y, rp2_z, 'Color',[0 0.45 0.74], 'LineWidth', 1.5);

plot3(rp2_x(1), rp2_y(1), rp2_z(1), '*b');

plot3(rp2_x(end), rp2_y(end), rp2_z(end), 'ob');

% plot final orbit

plot3(rp3_x, rp3_y, rp3_z, '-g', 'LineWidth', 1.5);

xlabel('X coordinate (ER)', 'FontSize', 12);

ylabel('Y coordinate (ER)', 'FontSize', 12);

zlabel('Z coordinate (ER)', 'FontSize', 12);

title('Initial, Transfer and Final Orbits', 'FontSize', 16);

axis equal;

grid on;

view(50, 20);

rotate3d on;

% create graphics disk file

print (gfilename1, '-dtiff');

saveas(h, gfilename2);


