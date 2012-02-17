fprintf('Animating %d timepoints... ', sp.Nt);

doQuiver = 0;

%% Normalize M
mz = squeeze(M(3,:,:,:)) / sp.Ms;
if doQuiver
    mx = squeeze(M(1,:,:,:)) / sp.Ms;
    my = squeeze(M(2,:,:,:)) / sp.Ms;
end

clf;
if doQuiver
    subplot(121);
end
ih = pcolor(mz(:,:,1));
shading interp;
%ih = imagesc(mz(:,:,1), [-1 1]);
%ih = imagesc(mz(:,:,1));
axis equal;
axis ij;
axis([1 sp.Nx 1 sp.Ny]);
ith = title('m_z(n = 0)');
colormap('hot');
% colormap('jet');
colorbar;

if doQuiver
    subplot(122);
    [X Y] = meshgrid([1:sp.Nx],[1:sp.Ny]); Z = zeros(size(X)); surf(X,Y,Z, 'facealpha',0.5, 'edgealpha',0); hold on;
    axis ij; axis equal; grid off; axis([1 sp.Nx 1 sp.Ny -1 1]);
    qh = quiver3(X,Y,Z, mx(:,:,1), my(:,:,1), mz(:,:,1)); %view(0,90);
    xlabel('x'); ylabel('y'); zlabel('z');
    view(35,30);
end


for i = 1:40:sp.Nt
    set(ih, 'cdata', double(mz(:,:,i)));
    set(ith, 'string', ['m_z(n = ', num2str(i), ')']);
    if doQuiver
        set(qh, 'udata', mx(:,:,i), 'vdata',my(:,:,i), 'wdata',mz(:,:,i));
    end
    % pause(0.01);
    drawnow;
end
fprintf('done\n');
