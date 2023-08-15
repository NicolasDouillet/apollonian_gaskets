function [] = apollonian_spherical_pseudo_gaskets(solid_id, nb_it, nb_samples)
% spherical_apollonian_pseudo_gaskets : function to compute and display the
% apollonian spherical pseudo gasket of one platonic solid (regular
% polyhedron).
%
% Author & support : nicolas.douillet (at) free.fr, 2023.


% Rotation matrices
Rmy = @(theta)[cos(theta)   0  -sin(theta);...
               0            1           0;...
               sin(theta)   0  cos(theta)];


Rmz = @(theta)[cos(theta) -sin(theta) 0;
               sin(theta)  cos(theta) 0;
               0           0          1];


switch solid_id
    
    case 1 % tetrahedron
                
        A = [2*sqrt(2)/3 0 -1/3]';
        B = [-sqrt(2)/3 sqrt(2)/sqrt(3) -1/3]';
        C = [0 0 1]';
        
        P = cat(2,A,B,C);
        
        R = apollonian_polygon_based_gasket(P, nb_it, nb_samples, false);        
        
        theta = 2*pi/3;
        u = A;
        
        Rmu = @(theta)[u(1)^2+cos(theta)*(1-u(1)^2) (1-cos(theta))*u(1)*u(2)-u(3)*sin(theta) (1-cos(theta))*u(1)*u(3)+u(2)*sin(theta);
                      (1-cos(theta))*u(1)*u(2)+u(3)*sin(theta) u(2)^2+cos(theta)*(1-u(2)^2) (1-cos(theta))*u(2)*u(3)-u(1)*sin(theta);
                      (1-cos(theta))*u(1)*u(3)-u(2)*sin(theta) (1-cos(theta))*u(2)*u(3)+u(1)*sin(theta) u(3)^2+cos(theta)*(1-u(3)^2)];
                
        R1 = zeros(size(R));
        R2 = zeros(size(R));
        R3 = zeros(size(R));
                
        for k = 1:size(R,3)
            
            R1(:,:,k) = Rmz(theta)*R(:,:,k);
            R2(:,:,k) = Rmz(2*theta)*R(:,:,k);
            R3(:,:,k) = Rmu(-theta)*R(:,:,k);
            
        end
                
        R = cat(3,R,R1,R2,R3);
        
    case 2 % cube / hexahedron
                
        P = [1  1  1  1;...
             1 -1 -1  1;...
             1  1 -1 -1];
        
        
        R = apollonian_polygon_based_gasket(P, nb_it, nb_samples, false);
        
        R1 = cat(3,R,cat(1,R(2,:,:),R(1,:,:),R(3,:,:)));
        R2 = cat(3,R1,cat(1,R(3,:,:),R(2,:,:),R(1,:,:)));
        R  = cat(3,R2,-R2);
        
    case 3 % octahedron                        
        
        P = eye(3);        
        R = apollonian_polygon_based_gasket(P, nb_it, nb_samples, false);        
        
        R1 = cat(3,R,cat(1,-R(1,:,:),R(2,:,:),R(3,:,:)));
        R1 = cat(3,R1,cat(1,R1(1,:,:),-R1(2,:,:),R1(3,:,:)));
        R  = cat(3,R1,-R1);
        
    case 4 % icosahedron
                
        [V,T] = platonic_solids(4,1,false,'default');                
        R = apollonian_polygon_based_gasket(V(T(3,:),:)',nb_it,nb_samples,false);
        R = cat(3,R,apollonian_polygon_based_gasket(V(T(11,:),:)',nb_it,nb_samples,false));
        
        % alpha = pi - asin(2/3); % icosahedron diedral angle
        
        n = size(R,3);
        R1 = R;
        R2 = R;
        R3 = R;
        R4 = R;
        
        for i = 1:n
            
            R1(:,:,i) = Rmz(0.4*pi)*R1(:,:,i);
            R2(:,:,i) = Rmz(0.8*pi)*R2(:,:,i);
            R3(:,:,i) = Rmz(1.2*pi)*R3(:,:,i);
            R4(:,:,i) = Rmz(1.6*pi)*R4(:,:,i);
            
        end
        
        R = cat(3,R,R1,R2,R3,R4);
        R = cat(3,R,-R);
        
    case 5 % dodecahedron
                
        [V,T] = platonic_solids(5,1,false,'default');                
        R = apollonian_polygon_based_gasket(V(T(7,:),:)',nb_it,nb_samples,false);
        
        alpha = pi - acos(-1/sqrt(5)); % dodecahedron diedral angle                
        
        n = size(R,3);
        R1 = zeros(size(R));
        R2 = zeros(size(R));
        R3 = zeros(size(R));
        R4 = zeros(size(R));
        R5 = zeros(size(R));
        
        for i = 1:n
            
            R1(:,:,i) = Rmz(0.2*pi)*Rmy(alpha)*R(:,:,i);
            
        end
        
        for i = 1:n
            
            R2(:,:,i) = Rmz(0.4*pi)*R(:,:,i);
            R3(:,:,i) = Rmz(0.8*pi)*R(:,:,i);
            R4(:,:,i) = Rmz(1.2*pi)*R(:,:,i);
            R5(:,:,i) = Rmz(1.6*pi)*R(:,:,i);
            
        end
        
        R = cat(3,R,R1,R2,R3,R4,R5);
        R = cat(3,R,-R);
        
        % otherwise
        
end


% projection on the sphere
N = sqrt(sum(R.^2,1));
R = R./N;

figure

for k = 1:size(R,3)
    
    line([R(1,:,k) R(1,1,k)],[R(2,:,k) R(2,1,k)],[R(3,:,k) R(3,1,k)],'Color',[0 0 1],'Linewidth',1), hold on;
    
end

axis equal, axis tight;
ax = gca;
ax.Clipping = 'off';


end % apollonian_spherical_pseudo_gaskets