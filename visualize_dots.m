%clear all;
%load ./matfiles/sim_10x10_dots_1e-08s_step_2e-12s_results.mat


fprintf('Animating %d timepoints\n', sp.Nt);

clf;
subplot(121); 
Ms = sp.Ms;
ih = imagesc(squeeze(M(3,:,:,1))/Ms, [-1 1]); 
axis equal; 
axis ij; 
axis([1 sp.Nx 1 sp.Ny]); 
ith = title('M_z(n = 0)');
colormap('hot'); 
%colorbar;

subplot(122);
Mx = squeeze(M(1,:,:,:));
My = squeeze(M(2,:,:,:));
Mz = squeeze(M(3,:,:,:));
[X Y] = meshgrid([1:sp.Nx],[1:sp.Ny]); Z = zeros(size(X)); surf(X,Y,Z, 'facealpha',0.5, 'edgealpha',0); hold on;
axis ij; axis equal; grid off; axis([1 sp.Nx 1 sp.Ny -1 1]);
qh = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), Mz(:,:,1)); %view(0,90);
xlabel('x'); ylabel('y'); zlabel('z');
view(0,90);



for i = 1:50:sp.Nt
    set(ih, 'cdata', double(squeeze(M(3,:,:,i))));
    set(ith, 'string', ['M_z(n = ', num2str(i), ')']);
    set(qh, 'udata', Mx(:,:,i), 'vdata',My(:,:,i), 'wdata',Mz(:,:,i));
    pause(0.01);
    drawnow;
end
