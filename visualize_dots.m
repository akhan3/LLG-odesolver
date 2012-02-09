%clear all;
%load ./matfiles/sim_10x10_dots_1e-08s_step_2e-12s_results.mat

%figure; 
clf;
colormap('jet'); 
%subplot(121); 
Ms = sp.Ms;
    ih = imagesc(squeeze(M(3,:,:,1)), [-Ms Ms]); axis equal; axis ij; axis([1 sp.Nx 1 sp.Ny]); 
    ith = title('M_z(t = 0)');
    %set(ih,'erasemode','background')
%     colorbar('westoutside'); 

%subplot(122); 
    %Mx = squeeze(M(1,:,:,:));
    %My = squeeze(M(2,:,:,:));
    %Mz = squeeze(M(3,:,:,:));
    %[X Y] = meshgrid([1:numdots_x],[1:numdots_y]); Z = zeros(size(X)); surf(X,Y,Z, 'facealpha',0.5, 'edgealpha',0); hold on; 
    %qth = title('Hext(t = 0)'); axis ij; axis equal; grid off; axis([1 numdots_x 1 numdots_y -1 1]); 
    %qh = quiver3(X,Y,Z, Mx(:,:,1), My(:,:,1), Mz(:,:,1)); %view(0,90);
    %xlabel('x'); ylabel('y'); zlabel('z'); 
    %%view(0,90);
    
%disp('paused!');
%keyboard;
%pause(1);
%disp('running...');


for i = 1:1:sp.Nt
%for i = 1:fieldlength
    set(ih, 'cdata', double(squeeze(M(3,:,:,i))));
    set(ith, 'string', ['M_z(t = ', num2str(t(i)), ')']);
    %set(qh, 'udata', Mx(:,:,i), 'vdata',My(:,:,i), 'wdata',Mz(:,:,i));
    %set(qth, 'string', ['Hext(t = ', num2str(t(i)), ')']);
    pause(0.01);
    drawnow;
end
