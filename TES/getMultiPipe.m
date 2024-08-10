function [VS_pipe,VS_boundaries,inds_pipe_boundary,inds_pipe_interior,boundary_perimiter_span] =  getMultiPipe(size_VS,x_center,y_center,ri,plot_bool)
if nargin == 4
    plot_bool = false;
end

Nx = size_VS(1);Ny = size_VS(2);

ii_center = floor(x_center)+1;
jj_center = floor(y_center)+1;
count = 0;
inds_pipe = [];
VS_pipe = false(Nx,Ny);
r_min = 1e6;
for dx_vx = round(-(ri+1)):round(ri+1)
    ii = ii_center + dx_vx;
    dx_cur = ii-0.5-x_center;
    for dy_vx = round(-(ri+1)):round(ri+1)
        jj = jj_center + dy_vx;
        dy_cur = jj-0.5-y_center;
        r_cur = (dx_cur^2 + dy_cur^2)^0.5;
        if r_cur <= ri
            count = count + 1;
            inds_pipe(count,:) = [ii,jj];
            VS_pipe(ii,jj) = true;
        end
        if r_cur < r_min
            r_min = r_cur;
        end
    end
end
N_pipe_cells = size(inds_pipe,1);
VS_boundaries = 6-countConnectivity(VS_pipe);
VS_boundaries(~VS_pipe) = 0;
N_boundary = sum(VS_boundaries(:)>0);

if N_pipe_cells == 1
    inds_pipe_boundary = inds_pipe;
    inds_pipe_interior = [];
    boundary_perimiter_span = 2*ri*pi;
else

    thetas_cell_boundary = zeros(N_boundary,1);
    inds_pipe_boundary = zeros(N_boundary,2);
    inds_pipe_interior = zeros(N_pipe_cells-N_boundary,2);
    count_boundary = 0;
    count_interior = 0;
    for nn = 1:N_pipe_cells
        ii = inds_pipe(nn,1);
        jj = inds_pipe(nn,2);
        if VS_boundaries(ii,jj)>0 % Then it is a boundary cell
            count_boundary = count_boundary+1;
            inds_pipe_boundary(count_boundary,:) = [ii,jj];
            thetas_cell_boundary(count_boundary) = atan2(jj-0.5-y_center,ii-0.5-x_center);
        else
            count_interior = count_interior+1;
            inds_pipe_interior(count_interior,:) = [ii,jj];
        end
    end
    [thetas_cell_boundary,sort_inds] = sort(thetas_cell_boundary);
    inds_pipe_boundary = inds_pipe_boundary(sort_inds,:);
    thetas_cell_boundary(end+1) = thetas_cell_boundary(1)+2*pi; % Make periodic
    diff_thetas = diff(thetas_cell_boundary); % Take difference
    thetas_cell_boundary(end) = [];
    thetas_cell_boundary = [thetas_cell_boundary(end)-2*pi;thetas_cell_boundary];
    diff_thetas = [diff(thetas_cell_boundary),diff_thetas];
    thetas_cell_boundary(1) = [];
    
    thetas_cell_edges = zeros(N_boundary,2);
    for nn = 1:N_boundary
        theta_cell = thetas_cell_boundary(nn);
        ii = inds_pipe_boundary(nn,1);
        jj = inds_pipe_boundary(nn,2);
    
        theta_b = theta_cell;
        theta_a = theta_b - 0.99*diff_thetas(nn,1);
        err = 1e6;
        while err>1e-6
            theta_mid = (theta_a + theta_b)/2;
            p_mid = [x_center,y_center] + [cos(theta_mid),sin(theta_mid)]*ri*0.9999999;
            ii_mid = floor(p_mid(1))+1;
            jj_mid = floor(p_mid(2))+1;
            [~,max_vec_ind] = max(abs([cos(theta_mid),sin(theta_mid)]));
            if max_vec_ind == 1 % x-component larger -> consider y boundaries
                vx_dif = abs(jj_mid-jj);
            else % consider x- boundaries
                vx_dif = abs(ii_mid-ii); 
            end
            if vx_dif > 0
                theta_a = theta_mid;
            else
                theta_b = theta_mid;
            end
            if vx_dif > 1
                err = 1e6;
            else
                err = min(abs(rem(p_mid,1)));
            end
        end
        thetas_cell_edges(nn,1) = theta_mid;
        theta_a = theta_cell;
        theta_b = theta_a + 0.99*diff_thetas(nn,2);
        err = 1e6;
        while err>1e-6
            theta_mid = (theta_a + theta_b)/2;
            p_mid = [x_center,y_center] + [cos(theta_mid),sin(theta_mid)]*ri*0.9999999; % Adjustment by factor of 0.9999999 prevents edge-cases on boundaries
            ii_mid = floor(p_mid(1))+1;
            jj_mid = floor(p_mid(2))+1;
            [~,max_vec_ind] = max(abs([cos(theta_mid),sin(theta_mid)]));
            if max_vec_ind == 1 % x-component larger -> consider y boundaries
                vx_dif = abs(jj_mid-jj);
            else % consider x- boundaries
                vx_dif = abs(ii_mid-ii);     
            end
            if vx_dif > 0
                theta_b = theta_mid;
            else
                theta_a = theta_mid;
            end
            if vx_dif > 1
                err = 1e6;
            else
                err = min(abs(rem(p_mid,1)));
            end
        end
        thetas_cell_edges(nn,2) = theta_mid;
    end
    theta_cell_edges_adj = zeros(N_boundary,2);
    for nn = 1:N_boundary
        nn_next = mod(nn,N_boundary)+1; % Circular indexing
        if abs(thetas_cell_edges(nn,2)-thetas_cell_edges(nn_next,1))>pi % They are separated by 2pi
            theta_l = min(thetas_cell_edges(nn,2),thetas_cell_edges(nn_next,1));
            theta_u = max(thetas_cell_edges(nn,2),thetas_cell_edges(nn_next,1));
            theta_shift = theta_l+2*pi;
            theta_mid = (theta_shift+theta_u)/2;
            if thetas_cell_edges(nn,2)<thetas_cell_edges(nn_next,1)
                thetas_cell_edges_adj(nn,2) = theta_mid-2*pi;
                thetas_cell_edges_adj(nn_next,1) = theta_mid;
            else
                thetas_cell_edges_adj(nn_next,1) = theta_mid-2*pi;
                thetas_cell_edges_adj(nn,2) = theta_mid;
            end
        else
            theta_mid = (thetas_cell_edges(nn,2) + thetas_cell_edges(nn_next,1))/2;
            thetas_cell_edges_adj(nn,2) = theta_mid;
            thetas_cell_edges_adj(nn_next,1) = theta_mid;
        end
    end
    boundary_perimiter_span = (thetas_cell_edges_adj(:,2)-thetas_cell_edges_adj(:,1))*ri;
    if abs(sum(boundary_perimiter_span)-2*pi*ri) > 1e-6
        error("Perimiter does not add up!")
    end
end 

if plot_bool
    f = figure;
    hold on;
    grid on;
    axis equal;
    for nn = 1:N_pipe_cells
        set(gca,'Color',[0.8,0.8,0.8])
        x_fill = inds_pipe(nn,1);
        y_fill = inds_pipe(nn,2);
        if ismember([x_fill,y_fill],inds_pipe_boundary,'rows')
            color = [0.9,0.9,0.9];
        else
            color = [1,1,1];

        end
        x_fill = [x_fill,x_fill,x_fill-1,x_fill-1];
        y_fill = [y_fill,y_fill-1,y_fill-1,y_fill];
        fill(x_fill,y_fill,color,'FaceAlpha',1,'LineWidth',2)  
    end
    if N_boundary ~= 1
        x_cell_edges = x_center + cos(thetas_cell_edges(:))*ri;
        y_cell_edges = y_center + sin(thetas_cell_edges(:))*ri;
        
        x_cell_edges_adj = x_center + cos(thetas_cell_edges_adj(:))*ri;
        y_cell_edges_adj = y_center + sin(thetas_cell_edges_adj(:))*ri;
        
        scatter(x_cell_edges_adj,y_cell_edges_adj,'o','MarkerEdgeColor','k','SizeData',20)
        if not(all(abs(x_cell_edges-x_cell_edges_adj)<1e-3) && all(abs(y_cell_edges-y_cell_edges_adj)<1e-3))
            scatter(x_cell_edges_adj,y_cell_edges_adj,'x','MarkerEdgeColor','k','MarkerFaceColor','k')
        end
    end
    all_thetas = linspace(-pi,pi,1000);
    all_x = x_center + cos(all_thetas)*ri;
    all_y = y_center + sin(all_thetas)*ri;
    plot(all_x,all_y,'Color',	"#0072BD",'LineWidth',1);
    xlim([x_center-round(ri)-1,x_center+round(ri)+1])
    ylim([y_center-round(ri)-1,y_center+round(ri)+1])
    xtick_bounds = xlim();
    xtick_bounds(1) = ceil(xtick_bounds(1));
    xtick_bounds(2) = floor(xtick_bounds(2));
    ytick_bounds = ylim();
    ytick_bounds(1) = ceil(ytick_bounds(1));
    ytick_bounds(2) = floor(ytick_bounds(2));
    xticks(xtick_bounds(1):1:xtick_bounds(2));
    yticks(ytick_bounds(1):1:ytick_bounds(2));
    set(gca,'GridColor',[0 0 0]);
    set(gca,'GridLineWidth',2)
   
end

end