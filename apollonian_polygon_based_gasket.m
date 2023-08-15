function [R, I, r] = apollonian_polygon_based_gasket(M, nb_it, nb_samples, option_display)
%% apollonian_polygon_based_gasket : function to compute and display the
% apollonian gasket of any 3D triangle or 3D regular polygon.
%
% Author & support : nicolas.douillet (at) free.fr, 2023.
%
%
% Syntax
%
% apollonian_polygon_based_gasket(P, nb_it, nb_samples);
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, option_display);
% [R, I, r] = apollonian_polygon_based_gasket(P, nb_it, nb_samples, option_display);
%
%
% Description
%
% apollonian_polygon_based_gasket(P, nb_it, nb_samples) compute and display the apollonian gasket of the polygon P
% at iteration #nb_it and using nb_samples for each circle.
%
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, option_display) displays the result (polygon + circles) when
% option_display is set to logical *true/1 (default), and doesn't when it is set to logical false/0.
%
% [R, I, r] = apollonian_polygon_based_gasket(P, nb_it, nb_samples, option_display)
% save the result (coordinates of every circles) in the array R.
%
%
% Input arguments
%
%       [-Mx-]
% - M = [-My-] : real matrix double. size(M,1) = 3. size(M,2) is the number of vertices of your polygon. 
%       [-Mz-]
%
% - nb_it : integer scalar, the number of iterations to perform. 
%
% - nb_samples : integer scalar, the number samples for each circle.
%
% - option_display : either logical, *true/false or numeric *1/0. 
%
%
% Output arguments
%
% - R : real matrix double. The apollonian circles coordinates. size(R) = [3, nb_samples, numel(r)], where
%       numel(r) is the number of resulting circles.
%
%       [-Ix-]
% - I = [-Iy-] : real matrix double. Coordinates of the sampled apollonian circles. size(R) = [3, numel(r)].
%       [-Iz-]
%
% r : real vector double. The apollonian circles radii / radius vector. 
%
%
% Example #1
% Equilateral triangle of the 3D space
% nb_it = 3;
% nb_samples = 64;
% P = eye(3);
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, true);
% view(70,21);
%
% Example #2
% Random triangle of the 3D space
% P = rand(3,3);
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, true);
%
% Example #3
% Square/rhombus of the 3D space
% V1 = [1 1 0]';
% V2 = [0 0 sqrt(2)]';
% V3 = [-1 -1 0]';
% V4 = [0 0 -sqrt(2)]';
% P = cat(2,V1,V2,V3,V4);
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, true);
% view(-132,-5);
%
% Example #4
% Pentagon of the 2D space
% P = [cos(0.4*pi) cos(0.8*pi) cos(1.2*pi) cos(1.6*pi) cos(2*pi);
%      sin(0.4*pi) sin(0.8*pi) sin(1.2*pi) sin(1.6*pi) sin(2*pi);
%      0           0           0           0           0        ];
%
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, true);
%
% Example #5
% Hexagon of the 3D space
% V1 = [-sqrt(2)/3 sqrt(6)/3 -1/3]';
% V3 = [-sqrt(2)/3 -sqrt(6)/3 -1/3]';
% V5 = [0 0 1]';
% V2 = (2/3)*(V1+V3) - (1/3)*V5;
% V4 = (2/3)*(V3+V5) - (1/3)*V1;
% V6 = (2/3)*(V1+V5) - (1/3)*V3;
% P = cat(2,V1,V2,V3,V4,V5,V6);
% apollonian_polygon_based_gasket(P, nb_it, nb_samples, true);


%% Body
nb_vertex = size(M,2);
nb_edges = nb_vertex;

uX = M - circshift(M,1,2); 

% Plane normal vector
n = cross(uX(:,1),uX(:,2));
n = n/norm(n);


% Compute first circle
if size(M,2) > 3 % polygon
    
    [R, I, r] = regular_3D_polygon_incircle(M',nb_samples,false);
    
elseif size(M,2) == 3 % triangle
    
    [R, I, r] = triangle_incircle(M(:,1),M(:,2),M(:,3),nb_samples,false);
    
    % else
    %
    %     error('Input polygon must have at least 3 vertices.');
    
end

R = permute(R,[2,1,3]);
I = I';


% Radial unitary vectors
rX = (M-I)./sqrt(sum((M-I).^2,1));


% Basic circle template
C_t = cat(1,cos(linspace(0,2*pi-2*pi/nb_samples,nb_samples)),...
            sin(linspace(0,2*pi-2*pi/nb_samples,nb_samples)),...
            zeros(1,nb_samples));

        
% Rotation in the polygon plane
% Vector u to rotate around
k = [0 0 1]';
u = cross(k,n);
alpha = atan2(norm(u),dot(k,n));
u = u/norm(u); % if norm(u) > 1e4*eps

% 3D rotation matrix around u vector
Rm = @(delta)[u(1,1)^2+cos(delta)*(1-u(1,1)^2) (1-cos(delta))*u(1,1)*u(2,1)-u(3,1)*sin(delta) (1-cos(delta))*u(1,1)*u(3,1)+u(2,1)*sin(delta);
              (1-cos(delta))*u(1,1)*u(2,1)+u(3,1)*sin(delta) u(2,1)^2+cos(delta)*(1-u(2,1)^2) (1-cos(delta))*u(2,1)*u(3,1)-u(1,1)*sin(delta);
              (1-cos(delta))*u(1,1)*u(3,1)-u(2,1)*sin(delta) (1-cos(delta))*u(2,1)*u(3,1)+u(1,1)*sin(delta) u(3,1)^2+cos(delta)*(1-u(3,1)^2)];

          
if alpha > 1e3*eps

    C_t = Rm(alpha) * C_t;
    
end


% Containers initialization
corner_idx      = {[]};
edge_idx        = {[]};
fathers_idx     = {[]}; % original / root circle has no father
sons_idx        = {[]};
soddys_done_idx = [];
      

for k = 1:nb_it
    
    
    circle_id = numel(r);
    
    
    for n = 1:circle_id                                
        
        
        if isempty(cell2mat(fathers_idx(n,1))) && isempty(cell2mat(sons_idx(n,1))) % Root circle first at iteration n = 1 only
            
            
            for corner_id = 1:nb_vertex % first nb_edges diagonal incircles
                
                [R, I, r] = compute_diagonal_incircle(R, I, r, M, rX, uX, circle_id, corner_id, nb_edges, nb_samples);
                
                idx = numel(r);
                
                fathers_idx(idx,1) = {1};
                sons_idx(idx,1)    = {[]};
                corner_idx(idx,1)  = {corner_id};
                edge_idx(idx,1)    = {[]};
                
            end
            
            sons_idx(1,1)   = {2:1+nb_vertex};
            
            
            % First Sangaku circles
            for circle_id = 2:1+nb_vertex                           
                
                for edge_id = [circle_id-1, 1+mod(circle_id-1,nb_edges)]         
                
                    [R, I, r] = compute_sangaku(R, I, r, M, C_t, uX, 1, circle_id, edge_id);
                    
                    idx = numel(r);
                    
                    fathers_idx(idx,1) = {[1,circle_id]};
                    sons_idx(idx,1)    = {[]};
                    corner_idx(idx,1)  = {[]};
                    edge_idx(idx,1)    = {edge_id};
                    sons_idx(circle_id,1)  = {[cell2mat(sons_idx(circle_id,1)),idx]};
                    
                end
                
            end
                        
            sons_idx(1,1)   = {[sons_idx{1}(1),2+nb_vertex:1+3*nb_vertex]};

            
            % First Soddy circles
            for circle_id1 = 2:1+nb_vertex                             
                
                for circle_id2 = 2*(circle_id1-1)+nb_vertex:2*(circle_id1-1)+nb_vertex+1
                    
                    [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, 1, circle_id1, circle_id2);
                    
                    idx = numel(r);
                    
                    fathers_idx(idx,1)          = {[1,circle_id1,circle_id2]};
                    sons_idx(idx,1)             = {[]};
                    corner_idx(idx,1)           = {[]};
                    edge_idx(idx,1)             = {[]};
                    sons_idx(circle_id1,1)      = {[cell2mat(sons_idx(circle_id1,1)),idx]};
                    sons_idx(circle_id2,1)      = {[cell2mat(sons_idx(circle_id2,1)),idx]};
                    
                end
                
            end
            
            sons_idx(1,1) = {[cell2mat(sons_idx(1,1)),2*+3*nb_vertex:1*+5*nb_vertex]};
            
        
        elseif numel(cell2mat(fathers_idx(n,1))) == 1 % diagonal circle                                                   
                        
            
            % Duplicata prevention            
            sons_ids = cell2mat(sons_idx(n,1));
            
            if ~isempty(sons_ids)
            
                sons_corner_ids = [];
                
                for j = sons_ids
                    
                    sons_corner_ids = cat(2,sons_corner_ids,cell2mat(corner_idx(j,1)));
                    
                end                               
                
            end
                                    
            
            if isempty(sons_ids) || (~isempty(sons_ids) && isempty(sons_corner_ids))
                
               [R, I, r] = compute_diagonal_incircle(R, I, r, M, rX, uX, n, corner_idx{n}(1), nb_edges, nb_samples);
               
               idx = numel(r);
               
               fathers_idx(idx,1) = {n};
               sons_idx(idx,1)    = {[]};
               corner_idx(idx,1)  = {corner_idx{n}(1)};
               edge_idx(idx,1)    = {[]};
               sons_idx(n,1)      = {[cell2mat(sons_idx(n,1)),idx]};
                
                
               % Linked Sangaku circles
               circle_id = numel(r);
               
               for edge_id = [cell2mat(corner_idx(circle_id,1)), 1+mod(cell2mat(corner_idx(circle_id,1)),nb_edges)]                   
                   
                   [R, I, r] = compute_sangaku(R, I, r, M, C_t, uX, fathers_idx{circle_id}(1), circle_id, edge_id);
                   
                   sgk_idx = numel(r);
                   
                   fathers_idx(sgk_idx,1) = {[cell2mat(fathers_idx(circle_id,1)),circle_id]};
                   sons_idx(sgk_idx,1)    = {[]};
                   corner_idx(sgk_idx,1)  = {[]};
                   edge_idx(sgk_idx,1)    = {edge_id};
                   sons_idx(n,1)          = {[cell2mat(sons_idx(n,1)),sgk_idx]};
                   
               end
               
                              
               % Two derived Soddy circles
               sgk_idx1 = numel(r) - 1;
               sgk_idx2 = numel(r);
               
               [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, fathers_idx{circle_id}(1), circle_id, sgk_idx1);
               [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, fathers_idx{circle_id}(1), circle_id, sgk_idx2);
               
               sd_idx1 = numel(r) - 1;
               sd_idx2 = numel(r);
               
               fathers_idx(sd_idx1,1) = {[fathers_idx{circle_id}(1),circle_id,sgk_idx1]};
               sons_idx(sd_idx1,1)    = {[]};
               corner_idx(sd_idx1,1)  = {[]};
               edge_idx(sd_idx1,1)    = {[]};
               
               fathers_idx(sd_idx2,1) = {[fathers_idx{circle_id}(1),circle_id,sgk_idx2]};
               sons_idx(sd_idx2,1)    = {[]};
               corner_idx(sd_idx2,1)  = {[]};
               edge_idx(sd_idx2,1)    = {[]};
               
               sons_idx(circle_id,1)                 = {[cell2mat(sons_idx(circle_id,1)),sd_idx1,sd_idx2]};
               sons_idx(fathers_idx{circle_id}(1),1) = {[cell2mat(sons_idx(fathers_idx{circle_id}(1),1)),sd_idx1,sd_idx2]};
               sons_idx(sgk_idx1,1)                  = {[cell2mat(sons_idx(sgk_idx1,1)),sd_idx1]};
               sons_idx(sgk_idx2,1)                  = {[cell2mat(sons_idx(sgk_idx2,1)),sd_idx2]};
               
                              
            end                                    
        
        
        elseif numel(cell2mat(fathers_idx(n,1))) == 2 % two fathers exactly (Sangaku circles)            
                        
            
            % Duplicata prevention
            sons_ids = cell2mat(sons_idx(n,1));
            
            if ~isempty(sons_ids)
                
                sons_edge_ids = [];
                
                for j = sons_ids
                    
                    sons_edge_ids = cat(2,sons_edge_ids,cell2mat(edge_idx(j,1)));
                    
                end
                
            end
            
            
            if isempty(sons_ids) || (~isempty(sons_ids) && isempty(sons_edge_ids))
                
                [R, I, r] = compute_sangaku(R, I, r, M, C_t, uX, fathers_idx{n}(1), n, edge_idx{n}(1));
                [R, I, r] = compute_sangaku(R, I, r, M, C_t, uX, fathers_idx{n}(2), n, edge_idx{n}(1));
                
                sgk_idx1 = numel(r) - 1;
                sgk_idx2 = numel(r);
                
                fathers_idx(sgk_idx1,1) = {[fathers_idx{n}(1),n]};
                sons_idx(sgk_idx1,1)    = {[]};
                corner_idx(sgk_idx1,1)  = {[]};
                edge_idx(sgk_idx1,1)    = {edge_idx{n}(1)};
                
                fathers_idx(sgk_idx2,1) = {[fathers_idx{n}(2),n]};
                sons_idx(sgk_idx2,1)    = {[]};
                corner_idx(sgk_idx2,1)  = {[]};
                edge_idx(sgk_idx2,1)    = {edge_idx{n}(1)};
                
                sons_idx(fathers_idx{n}(1),1) = {[cell2mat(sons_idx(fathers_idx{n}(1),1)),sgk_idx1]};
                sons_idx(fathers_idx{n}(2),1) = {[cell2mat(sons_idx(fathers_idx{n}(2),1)),sgk_idx2]};
                sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),sgk_idx1,sgk_idx2]};
                
                
                % Derived Soddy circles
                [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, n, fathers_idx{n}(1), sgk_idx1);
                [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, n, fathers_idx{n}(2), sgk_idx2);
                
                sd_idx1 = numel(r) - 1;
                sd_idx2 = numel(r);
                
                fathers_idx(sd_idx1,1) = {[n,fathers_idx{n}(1),sgk_idx1]};
                sons_idx(sd_idx1,1)    = {[]};
                corner_idx(sd_idx1,1)  = {[]};
                edge_idx(sd_idx1,1)    = {[]};
                
                sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),sd_idx1]};
                sons_idx(fathers_idx{n}(1),1) = {[cell2mat(sons_idx(fathers_idx{n}(1),1)),sd_idx1]};
                sons_idx(sgk_idx1,1)          = {[cell2mat(sons_idx(sgk_idx1,1)),sd_idx1]};
                
                
                fathers_idx(sd_idx2,1) = {[n,fathers_idx{n}(2),sgk_idx2]};
                sons_idx(sd_idx2,1)    = {[]};
                corner_idx(sd_idx2,1)  = {[]};
                edge_idx(sd_idx2,1)    = {[]};
                
                sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),sd_idx2]};
                sons_idx(fathers_idx{n}(2),1) = {[cell2mat(sons_idx(fathers_idx{n}(2),1)),sd_idx2]};
                sons_idx(sgk_idx2,1)          = {[cell2mat(sons_idx(sgk_idx2,1)),sd_idx2]};
                
                
            end
                       
            
        elseif numel(cell2mat(fathers_idx(n,1))) == 3 && isempty(find(soddys_done_idx == n, 1)) % three fathers exactly : Soddy circles
                        
            
            [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, fathers_idx{n}(1), fathers_idx{n}(2), n);
            [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, fathers_idx{n}(1), n, fathers_idx{n}(3));
            [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, n, fathers_idx{n}(2), fathers_idx{n}(3));
            
            new_idx = numel(r);
            
            
            % father #1
            fathers_idx(new_idx-2,1) = {[fathers_idx{n}(1),fathers_idx{n}(2),n]};
            sons_idx(new_idx-2,1)    = {[]};
            corner_idx(new_idx-2,1)  = {[]};
            edge_idx(new_idx-2,1)    = {[]};
            
            sons_idx(fathers_idx{n}(1),1) = {[cell2mat(sons_idx(fathers_idx{n}(1),1)),new_idx-2]};
            sons_idx(fathers_idx{n}(2),1) = {[cell2mat(sons_idx(fathers_idx{n}(2),1)),new_idx-2]};
            sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),new_idx-2]};
            
            
            % father #2
            fathers_idx(new_idx-1,1) = {[fathers_idx{n}(1),n,fathers_idx{n}(3)]};
            sons_idx(new_idx-1,1)    = {[]};
            corner_idx(new_idx-1,1)  = {[]};
            edge_idx(new_idx-1,1)    = {[]};
            
            sons_idx(fathers_idx{n}(1),1) = {[cell2mat(sons_idx(fathers_idx{n}(1),1)),new_idx-1]};
            sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),new_idx-1]};
            sons_idx(fathers_idx{n}(3),1) = {[cell2mat(sons_idx(fathers_idx{n}(3),1)),new_idx-1]};
            
            
            % father #3
            fathers_idx(new_idx,1) = {[n,fathers_idx{n}(2),fathers_idx{n}(3)]};
            sons_idx(new_idx,1)    = {[]};
            corner_idx(new_idx,1)  = {[]};
            edge_idx(new_idx,1)    = {[]};
            
            sons_idx(n,1)                 = {[cell2mat(sons_idx(n,1)),new_idx]};
            sons_idx(fathers_idx{n}(2),1) = {[cell2mat(sons_idx(fathers_idx{n}(2),1)),new_idx]};
            sons_idx(fathers_idx{n}(3),1) = {[cell2mat(sons_idx(fathers_idx{n}(3),1)),new_idx]};
            
            
            % Store id for duplicata prevention
            soddys_done_idx = [soddys_done_idx, n];

        end
        
        
    end
    
    
end

    
%% Display
if option_display
    
    figure        
    line([M(1,:) M(1,1)],[M(2,:) M(2,1)],[M(3,:) M(3,1)],'Color',[1 0 0],'Linewidth',2), hold on;    
    
    for k = 1:size(R,3)
        
        line([R(1,:,k) R(1,1,k)],[R(2,:,k) R(2,1,k)],[R(3,:,k) R(3,1,k)],'Color',[0 0 1],'Linewidth',1), hold on;
        
    end        
    
    view(3);    
    axis equal, axis tight;    
    ax = gca;
    ax.Clipping = 'off';
    
end


end % apollonian_polygon_based_gasket


%% compute_diagonal_incircle subfunction
function [R, I, r] = compute_diagonal_incircle(R, I, r, M, rX, uX, circle_id, corner_id, nb_edges, nb_samples)


In = I(:,circle_id);
rn = r(circle_id);

uX1 = uX(:,corner_id);
uX2 = uX(:,1+mod(corner_id,nb_edges));
rX = rX(:,corner_id);

Mid = M(:,corner_id);
xK = In + rn*rX;

% New vertices
U = line_plane_intersection(uX2,Mid,rX,xK);
V = line_plane_intersection(uX1,Mid,rX,xK);

% Circle inscribed in the new triangle
[cM,Im,Rm] = triangle_incircle(Mid,U,V,nb_samples,false);

I = cat(2,I,Im');
r = cat(2,r,Rm);
R = cat(3,R,cM');


end % compute_diagonal_incircle


%% compute_sangaku subfunction
function [R, I, r] = compute_sangaku(R, I, r, M, C_t, uX, father1_idx, father2_idx, edge_id)


Mid = M(:,edge_id);
uX = uX(:,edge_id);

% New radius
rSX = r(1,father2_idx) * r(1,father1_idx) / (sqrt(r(1,father2_idx)) + sqrt(r(1,father1_idx)))^2;

% New centres
[~,HMX] = point_to_line_distance(I(:,father1_idx)',uX',Mid');
[~,HIM] = point_to_line_distance(I(:,father2_idx)',uX',Mid');

hXI = (I(:,father1_idx)-HMX')/norm((I(:,father1_idx)-HMX'));
uIXI = (HIM'-HMX')/norm(HIM'-HMX');

% New circle coordinates
IsX = HMX' + 2*sqrt(r(1,father1_idx)*rSX)*uIXI + rSX*hXI;
RsX = IsX + rSX*C_t;

R = cat(3,R,RsX);
I = cat(2,I,IsX);
r = cat(2,r,rSX);


end % compute_sangaku


%% compute_soddy_circle subfunction
function [R, I, r] = compute_soddy_circle(R, I, r, Rm, alpha, C_t, circle_id1, circle_id2, circle_id3)


X1 = I(:,circle_id1);
X2 = I(:,circle_id2);
X3 = I(:,circle_id3);

r1 = r(circle_id1);
r2 = r(circle_id2);
r3 = r(circle_id3);

gamma1 = 1/r1;
gamma2 = 1/r2;
gamma3 = 1/r3;

gamma_S = gamma1 + gamma2 + gamma3 + 2*sqrt(gamma1*gamma2 + gamma2*gamma3 + gamma3*gamma1);
rS = 1/gamma_S;

% Rotation of I1, I2, and I3 in (Oxy) plane
G = mean([I(:,circle_id1),I(:,circle_id2),I(:,circle_id3)],2);

I1 = Rm(-alpha)*(X1-G);
I2 = Rm(-alpha)*(X2-G);
I3 = Rm(-alpha)*(X3-G);

z1 = I1(1,1) + 1i*I1(2,1);
z2 = I2(1,1) + 1i*I2(2,1);
z3 = I3(1,1) + 1i*I3(2,1);

% 2D solve in this plane
zS1 = (gamma1*z1 + gamma2*z2 + gamma3*z3 + 2*sqrt(gamma1*gamma2*z1*z2 + gamma2*gamma3*z2*z3 + gamma3*gamma1*z3*z1)) / gamma_S;
zS2 = (gamma1*z1 + gamma2*z2 + gamma3*z3 - 2*sqrt(gamma1*gamma2*z1*z2 + gamma2*gamma3*z2*z3 + gamma3*gamma1*z3*z1)) / gamma_S;

Is1 = cat(1,real(zS1),imag(zS1),0);
Is2 = cat(1,real(zS2),imag(zS2),0);
    
Is1 = G + Rm(alpha)*Is1;
Is2 = G + Rm(alpha)*Is2;


if abs(zS1) > abs(zS2)
    
    Is = Is1;
    
else
    
    Is = Is2;
    
end

Rs = Is + rS*C_t;

r = cat(2,r,rS);
I = cat(2,I,Is);
R = cat(3,R,Rs);


end % compute_soddy_circle