%% animation for group 3d gradient

filename = 'Group_Animated.gif';

	n = 0;
    % 3d scatter plot, coloured by yeo
		h = figure('units','centimeter','outerposition',[0 0 80 80]);
		a(1) = axes('position', [0.05 0.3 0.4 0.4]);
  	scatter3(group_grad(:,1), group_grad(:,2), group_grad(:,3), 250, yeo_parc, 'filled');
		grid off
		xlabel('G1','FontSize', 24);
		ylabel('G2','FontSize', 24);
		zlabel('G3','FontSize', 24);
		set(a(1), 'Color', [1,1,1]);
		set(gcf, 'Color', [1,1,1]);
		set(gcf, 'InvertHardCopy', 'off');
		set(gca, 'Visible','off');
    colormap(a(1), cmap_yeo)
    view(n,20)
    gif(filename, 'DelayTime', 1/6)
	for n = 1:3:360
		view(n,20)
		gif
	end


%% animation for age change
age_decile = (age>=18) + (age>=28) + (age>=38) + (age>=48) + (age>=58) + (age>=68) + (age>=78) + (age>=88);
age_decile_title = {'18-27' '28-37' '38-47' '48-57' '58-67' '68-77' '78-87'};

			filename = '2D_Animated.gif';
	    % plot change in 2d space

	    combos = combnk(1:3,2);
			for ii = 1:max(age_decile)
				% get mean CoG
				decile_yeo_cog(ii,:,:) = mean(yeo_cog(age_decile==ii,:,:),1);
				h = figure('units','centimeter','outerposition',[0 0 30 60]);
	    	for c = 1:length(combos)
					decile_grad_x = mean(G1(:,combos(c,1),age_decile==ii),3);
					decile_grad_y = mean(G1(:,combos(c,2),age_decile==ii),3);
	        a(c+1) = axes('position', [0.3 1-(0.3*c) 0.5 0.26]);
	        scatter( decile_grad_x, decile_grad_y, 100, yeo_parc, 'filled', 'MarkerFaceAlpha', 0.2);
					hold on;
	        xlim([-0.1 0.15])
	        ylim([-0.1 0.1])
	        colormap(a(c+1), cmap_yeo)
	        scatter(squeeze(decile_yeo_cog(ii,combos(c,1),:)), squeeze(decile_yeo_cog(ii,combos(c,2),:)), 250, cmap_yeo(2:end,:), 'filled');
					set(gca, 'Visible','off');
					set(gcf, 'Color', [1,1,1]);
					set(gcf, 'InvertHardCopy', 'off');
		  end
			hold off
			%suptitle(age_decile_title{ii});
			axes( 'Position', [0, 0.9, 0.3, 0.3] ) ;
			set( gca, 'Color', 'None', 'XColor', 'White', 'YColor', 'White' ) ;
			text(0.2, 0, age_decile_title{ii}', 'FontSize', 24', 'FontWeight', 'Bold') ;
			frame = getframe(h);
			im = frame2im(frame);
			[imind,cm] = rgb2ind(im,256);

			if ii == 1
				 imwrite(imind,cm,filename,'gif', 'Loopcount',inf,'DelayTime', 2);
			else
				imwrite(imind,cm,filename,'gif','WriteMode','append','DelayTime', 2);
			end

		end
