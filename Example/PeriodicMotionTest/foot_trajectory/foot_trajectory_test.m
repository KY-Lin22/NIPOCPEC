% phi \in (-pi, pi) varies in a clock-wise manner


%% case 1 (stance): phi from pi to 0, then p from 1 to -1
clear all
clc

stages = 100;
phi_st = pi : -(pi/stages) : 0;
p_st = 2*phi_st/pi - 1;

T = 0.5;
beta = 1/2;
vx = 0.5;
ax = vx * T * beta;

rx_st = zeros(1, stages + 1);
for n = 1 : stages + 1
    rx_st(n) = ax * (6 * p_st(n)^5 - 15 * p_st(n)^4 + 10 * p_st(n)^3 - 0.5);
end
rz_st = zeros(1, stages + 1);

% prepare for animation
phi = phi_st;
rx = rx_st;
rz = rz_st;

%% case 2 (swing): phi from 0 to -pi, then p from -1 to 1
clear all
clc

stages = 100;
phi_sw = 0 : -(pi/stages) : -pi;
p_sw = -2*phi_sw/pi - 1;

T = 0.5;
beta = 1/2;
vx = 0.5;
ax = vx * T * beta;

rx_sw = zeros(1, stages + 1);
for n = 1 : stages + 1
    rx_sw(n) = ax * (6 * p_sw(n)^5 - 15 * p_sw(n)^4 + 10 * p_sw(n)^3 - 0.5);
end
h = 0.08;
rz_sw = zeros(1, stages + 1);
for n = 1 : stages + 1
    rz_sw(n) = h * (-64 * p_sw(n)^6 + 192 * p_sw(n)^5 - 192*p_sw(n)^4 + 64*p_sw(n)^3);
end

% prepare for animation
phi = phi_sw;
rx = rx_sw;
rz = rz_sw;

%%

traj_X = cell(1, stages + 1);
traj_Y = cell(1, stages + 1);
for n = 1 : stages + 1
    traj_X{1, n} = [rx(1, 1 : n), repmat(rx(1, n), 1, stages + 1 - n)];
    traj_Y{1, n} = [rz(1, 1 : n), repmat(rz(1, n), 1, stages + 1 - n)];
end
% get figure size
figure(100)
pos = get(gcf, 'Position');
width = pos(3);
height = pos(4);
% define movie record
mov = zeros(height*3/2, width*3/2, 1, stages + 1, 'uint8');

% pre allocate
traj = plot(traj_X{1, 1}, traj_Y{1, 1}, '-', 'Color', [1 0 0], 'MarkerSize', 1, 'LineWidth', 1);
axisLimit_X = [min(rx) - 0.5; max(rx) + 0.5];
axisLimit_Y = [min(rz) - 0.5; max(rz) + 0.5]; 
axis([axisLimit_X; axisLimit_Y]);
phiPrint = title(sprintf('phi: %0.2f', phi(1)));
xlabel('x')
ylabel('z')

for n = 1 : stages + 1
    % update XData and YData   
    set(traj, 'XData', traj_X{1, n}, 'YData', traj_Y{1, n});
    set(phiPrint, 'String', sprintf('phi: %0.2f sec', phi(n)));
    % get frame as an image
    f = getframe(gcf);
    % Create a colormap for the first frame. for the rest of the frames, use the same colormap
    if n == 1
        [mov(:,:,1,n), map] = rgb2ind(f.cdata, 256, 'nodither');
    else
        mov(:,:,1,n) = rgb2ind(f.cdata, map, 'nodither');
    end
end
% create an animated GIF
imwrite(mov, map, 'gait.gif', 'DelayTime', 0, 'LoopCount', inf)