function plot_nodes_direction(nodefile, facefile)

% Figure
hf = figure(1);
set(hf, 'Position', [100 100 600 600]);

% Plot nodes
nodes = load(nodefile);

plot3(nodes(:,1), nodes(:,2), nodes(:,3), 'k.');

grid;

% Plot faces
faces = load(facefile);

xnodes = nodes(:,1);
ynodes = nodes(:,2);
znodes = nodes(:,3);

xfaces = xnodes(faces+1);
yfaces = ynodes(faces+1);
zfaces = znodes(faces+1);

for i = 1:size(xfaces,1)
    total = 0;
    for k = 1:4
        kx = mod(k,4)+1;
        lat1 = asin(znodes(faces(i,k)+1));
        lat2 = asin(znodes(faces(i,kx)+1));

        lon1 = atan2(ynodes(faces(i,k)+1), xnodes(faces(i,k)+1));
        lon2 = atan2(ynodes(faces(i,kx)+1), xnodes(faces(i,kx)+1));
        
        if ((lon1 > pi/2) && (lon2 < -pi/2))
            lon2 = lon2 + 2*pi;
        end
        if ((lon1 < -pi/2) && (lon2 > pi/2))
            lon1 = lon1 + 2*pi;
        end
        
        total = total + (lon2 - lon1) * (lat2 + lat1);
    end

    if (total >= 0)
        patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'r');
    else
        patch(xfaces(i,:), yfaces(i,:), zfaces(i,:),'w');
    end
end

for i = 1:size(xfaces,1)
    total = 0;
    for k = 1:4
        kx = mod(k,4)+1;
        lat1 = asin(znodes(faces(i,k)+1));
        lat2 = asin(znodes(faces(i,kx)+1));

        lon1 = atan2(ynodes(faces(i,k)+1), xnodes(faces(i,k)+1));
        lon2 = atan2(ynodes(faces(i,kx)+1), xnodes(faces(i,kx)+1));

        if ((lon1 > pi/2) && (lon2 < -pi/2))
            lon2 = lon2 + 2*pi;
        end
        if ((lon1 < -pi/2) && (lon2 > pi/2))
            lon1 = lon1 + 2*pi;
        end

        total = total + (lon2 - lon1) * (lat2 + lat1);
    end
    
    if (total < 0)
        arrow3d([xfaces(i,1) yfaces(i,1) zfaces(i,1)], [xfaces(i,2) yfaces(i,2) zfaces(i,2)]);
        arrow3d([xfaces(i,2) yfaces(i,2) zfaces(i,2)], [xfaces(i,3) yfaces(i,3) zfaces(i,3)]);
        arrow3d([xfaces(i,3) yfaces(i,3) zfaces(i,3)], [xfaces(i,4) yfaces(i,4) zfaces(i,4)]);
        arrow3d([xfaces(i,4) yfaces(i,4) zfaces(i,4)], [xfaces(i,1) yfaces(i,1) zfaces(i,1)]);
        disp(sprintf('Face: %i', i));
        break;
    end
end

set(gca, 'FontSize', 14);