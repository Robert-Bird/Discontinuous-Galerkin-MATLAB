function [vert,conn,tria] = triangle_mesh(node,edge,hmax)
P = node(:,1:2);                                                           % point list (original nodes kept as rows 1:n0)
C = zeros(0,2);                                                            % boundary constraint segments, indices into P

% (1) Resample each boundary edge at ~hmax -------------------------------
for e = 1:size(edge,1)
    a = edge(e,1); b = edge(e,2);
    m = max(1,round(hypot(P(a,1)-P(b,1),P(a,2)-P(b,2))/hmax));             % number of sub-segments on this edge
    if m == 1
        C = [C; a b];                                                      %#ok<AGROW>
    else
        t     = (1:m-1).'/m;
        idx   = size(P,1) + (1:m-1).';                                     % indices of the new edge points
        P     = [P; (1-t).*P(a,:) + t.*P(b,:)];                            %#ok<AGROW>
        chain = [a; idx; b];                                               % a -> interior points -> b
        C     = [C; chain(1:end-1) chain(2:end)];                          %#ok<AGROW>
    end
end

% (2) Seed interior grid points at ~hmax, kept inside and off the boundary
DTb = delaunayTriangulation(P,C);                                          % boundary-only triangulation, used to test points
[gx,gy] = meshgrid(min(P(:,1)):hmax:max(P(:,1)), min(P(:,2)):hmax:max(P(:,2)));
gp   = [gx(:) gy(:)];
tID  = pointLocation(DTb,gp);                                              % which triangle each grid point falls in (NaN if outside hull)
inDT = isInterior(DTb);
keep = ~isnan(tID);
keep(keep) = inDT(tID(keep));                                              % inside the (possibly non-convex) domain
dmin = inf(size(gp,1),1);                                                  % distance from each grid point to the nearest boundary point
for j = 1:size(P,1)                                                        % (P is boundary points only at this stage)
    dmin = min(dmin,hypot(gp(:,1)-P(j,1),gp(:,2)-P(j,2)));
end
keep = keep & dmin > 0.5*hmax;                                             % keep grid points clear of the boundary
P = [P; gp(keep,:)];

% (3) Triangulate the full point set -------------------------------------
DT     = delaunayTriangulation(P,C);
inside = isInterior(DT);
tria   = DT.ConnectivityList(inside,:);
vert   = DT.Points;
conn   = freeBoundary(triangulation(tria,vert));

end
