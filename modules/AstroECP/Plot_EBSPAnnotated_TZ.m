function [ f1 ] = Plot_EBSPAnnotated_TZ( EBSD_Pattern_in,EBSD_Geometry,~,RotMat,Crystal_UCell,Family_List,s1, keV)

% PLOT_EBSPAnnotated_TZ
% Plot an annotated EBSP with band centre and edges
% Created by Tianbi Zhang on Feb 28 2022 based on original code
% Last modified on Mar 01 2022 by Tianbi Zhang
% Code is annotated by Tianbi Zhang on Mar 01 2022

%% Initialization
% Color sets
cset=[1,0,0;
    0,1,0;
    0,0,1;
    0,1,1;
    1,0,1;
    1,1,0];
cset=[cset;cset*0.4];

%U.K         =   Crystal_UCell.At*RotMat;
%U.Kstar     =   Crystal_UCell.Astar*RotMat;

% Calculate electron wavelength
% V = 20 * 1000; % electron voltage in V
V = keV; % electron voltage in V

me = 9.1 * 10^(-31); %electron mass in kg
e_charge = 1.6 * 10^(-19); % unit charge in C
c = 3 * 10^8; % Speed of light m/s
h = 6.626 * 10^(-34); % Planck constant in Js

% Incorporate relativistic correction
% based on: https://advanced-microscopy.utah.edu/education/electron-micro
lambda = ( h / sqrt(2 * me * e_charge * V) ) / (sqrt(1 + e_charge * V /2/me/c^2)); % in m

% The following two lines only apply to cubic materials!!!
% For non-cubic systems, please consider using a more generalized formula.
cubic_a = Crystal_UCell.a; % Extract lattice constant in ang from unit cell variable
cubic_a = cubic_a / 10^(10); % convert to m

%% Plot the original EBSP
f1(1)=imagesc(EBSD_Geometry.x_screen,EBSD_Geometry.y_screen,EBSD_Pattern_in,'Parent',s1);

axis(s1,'equal','xy');
hold(s1,'on')
%label the axes
xlabel(s1,'X / Z');
ylabel(s1,'Y / Z');

colormap('gray')

% % plot the PC
% f1(2) = scatter(0, 0, 150, 'wo', 'filled', 'MarkerEdgeColor', 'k','Parent', s1);
% f1(3) = scatter(0, 0, 150, 'kx', 'Parent', s1);



%% Annotate band centre and edges
% plot the bands from the rotated solution

if ~isempty(Family_List)
    for n=1:size(Family_List,1)
        Plane_HKLs=Family_List{n}(:,2:4); % reads indices of plane reflectors
        HKLs=Family_List{n}(:,5:7)*RotMat; % reads indices of plane normals
        
        for p=1:size(HKLs,1)
            nhat=HKLs(p,:); % read the plane normal and assign it to nhat
            nhat=nhat/norm(nhat); % normalize to a unit vector
            
            % calculate the d-spacing and Bragg angle
            d_hkl = cubic_a / norm(Plane_HKLs(p,:)); % N.B. This only works for cubic! More generalized formula needed for hex, etc...
            sin_theta = lambda / 2 / d_hkl; % Bragg angle
            
            % Below plot the centre line and edges based on the orientation
            % of the plane normal
            if abs(nhat(2)) > 0.5
                xlin=linspace(EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end),100); % vector that spans the entire x axis
                
                % Plot band centre utilizing the relationship 
                % nhat \cdot l = 0 and Z=1.
                %ylin=(-xlin.*nhat(1)-nhat(3))./nhat(2); % For band centre, utilize nhat \cdot l = 0 and z=1
                
%                 Choose only the part of Y within the EBSP frame
%                 [ix]=find(ylin > EBSD_Geometry.y_screen(1) & ylin < EBSD_Geometry.y_screen(end)); % 
%                 ylin2=ylin(ix);
%                 xlin2=xlin(ix);
                
                % Plot band edges utilizing the relationship
                % nhat \cdot l = sin(theta_(hkl)) where theta is the Bragg
                % angle.
                y_edge1_lin = (-xlin.*nhat(1)-nhat(3)+sin_theta)./nhat(2);
                y_edge2_lin = (-xlin.*nhat(1)-nhat(3)-sin_theta)./nhat(2);

                % Plot band centre (solid line) and edges (dashed line)
%                 if ~isempty(ix)
%                     f1(end+1)=plot(xlin,ylin,'color',cset(n,:),'LineWidth',2,'parent',s1);
                    f1(end+1)=plot(xlin,y_edge1_lin,'color',cset(n,:),'LineWidth',0.5,'parent',s1,'LineStyle','-');
                  % f1(end+1)=plot(xlin,y_edge2_lin,'color',cset(n,:),'LineWidth',1,'parent',s1,'LineStyle','-');

%                 end
            else
                % below is the same calculation for another case of nhat
                ylin=linspace(EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end),100);
                xlin=(-ylin.*nhat(2)-nhat(3))./nhat(1);
                
%                 [ix]=find(xlin > EBSD_Geometry.x_screen(1) & xlin < EBSD_Geometry.x_screen(end));
%                 ylin2=ylin(ix);
%                 xlin2=xlin(ix);
                
                x_edge1_lin = (-ylin.*nhat(2)-nhat(3)+sin_theta)./nhat(1);
                x_edge2_lin = (-ylin.*nhat(2)-nhat(3)-sin_theta)./nhat(1);

%                 if ~isempty(ix)
%                     f1(end+1)=plot(xlin2,ylin2,'color',cset(n,:),'LineWidth',2,'Parent',s1);
                    f1(end+1)=plot(x_edge1_lin,ylin,'color',cset(n,:),'LineWidth',0.5,'parent',s1,'LineStyle','-');
                   % f1(end+1)=plot(x_edge2_lin,ylin,'color',cset(n,:),'LineWidth',2,'parent',s1,'LineStyle','-');
%                 end
                
            end
            
        end
    end
end

%plot the detected bands
% Note that for ECCI_Si this will not do anything since input nhat_gnom = [].
% if ~isempty(nhat_gnom)
%     for n=1:size(nhat_gnom,1)
%         
%         nhat=nhat_gnom(n,:);
%         nhat=nhat/norm(nhat);
%         
%         if abs(nhat(2)) > 0.5
%             xlin=linspace(EBSD_Geometry.x_screen(1),EBSD_Geometry.x_screen(end),100);
%             ylin=(-xlin.*nhat(1)-nhat(3))./nhat(2);
%             
%             [ix]=find(ylin > EBSD_Geometry.y_screen(1) & ylin < EBSD_Geometry.y_screen(end));
%             ylin2=ylin(ix);
%             xlin2=xlin(ix);
%             
%             if ~isempty(ix)
%                 f1(end+1)=plot(xlin2,ylin2,'.','color','w','LineWidth',2,'Parent',s1);
%             end
%         else
%             
%             ylin=linspace(EBSD_Geometry.y_screen(1),EBSD_Geometry.y_screen(end),100);
%             xlin=(-ylin.*nhat(2)-nhat(3))./nhat(1);
%             
%             [ix]=find(xlin > EBSD_Geometry.x_screen(1) & xlin < EBSD_Geometry.x_screen(end));
%             ylin2=ylin(ix);
%             xlin2=xlin(ix);
%             
%             if ~isempty(ix)
%                 f1(end+1)=plot(xlin2,ylin2,'.','color','w','LineWidth',2,'Parent',s1);
%             end
%             
%         end
%         
%     end
% end


% Adjust X, Y limits in the end to match that of the original pattern
ylim(s1,[EBSD_Geometry.y_screen(1) EBSD_Geometry.y_screen(end)]);
xlim(s1,[EBSD_Geometry.x_screen(1) EBSD_Geometry.x_screen(end)]);
end

