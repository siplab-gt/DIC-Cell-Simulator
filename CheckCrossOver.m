function cross_itself = CheckCrossOver(x_temp)

% This function checks whether coordinates in x_temp crosses over itself.
% Inputs:
% x_temp            Coordinates organized as (x1,x2,...,xK,y1,y2,...,yK)
% Outputs:
% cross_itself      Boolean value that denotes whether there is cross-over.

d = length(x_temp);

cross_itself = false;
for p = 1:d/2
    % reshuffle points such that p-th point is the first coordinate
    x_temp(1:d/2) = [ x_temp(p:d/2); x_temp(1:p-1) ];
    x_temp(d/2+1:end) = [ x_temp([p:d/2]+d/2); x_temp([1:p-1]+d/2) ];
    % normalize by first coordinate
    x_temp(1:d/2) = x_temp(1:d/2) - x_temp(1);
    x_temp(d/2+1:end) = x_temp(d/2+1:end) - x_temp(d/2+1);
    % rotate such that coord1 and coord2 is horizontal
    rot_angle = -atan2( x_temp(d/2+2), x_temp(2) );
    rot_matrix = [cos(rot_angle),-sin(rot_angle); sin(rot_angle), cos(rot_angle)];
    x_rot = rot_matrix * [ x_temp(1:d/2), x_temp(d/2+1:end)]';
    x_rot = [ reshape(x_rot(1,:),[],1); reshape(x_rot(2,:),[],1) ];
    % rotate entire cell by angle from two coordinates
    coord1 = [ x_rot(1); x_rot(d/2+1) ];
    coord2 = [ x_rot(2); x_rot(d/2+2) ];
    % consider all other points
    for k = 3:d/2-1
        coord3 = [ x_rot(k); x_rot(d/2+k) ];
        coord4 = [ x_rot(k+1); x_rot(d/2+k+1) ];
        m1 = ( coord2(2) - coord1(2) ) / ( coord2(1) - coord1(1) );
        c1 = coord2(2) - m1 * coord2(1);
        m2 = ( coord4(2) - coord3(2) ) / ( coord4(1) - coord3(1) );
        c2 = coord4(2) - m2 * coord4(1);
        x_cross = ( c2 - c1 ) / (m1 - m2); % analytically determine crossing
        cross_itself = cross_itself | ...
                       ( (x_cross > coord1(1)) & (x_cross < coord2(1) ) & ...
                         ( (x_cross > coord3(1)) & (x_cross < coord4(1) ) | ...
                           (x_cross > coord4(1)) & (x_cross < coord3(1) ) ) );
    end
    if cross_itself == true, break; end;
end

end