function [seg, pbar, F, inner] = airfoil_discretize(datfile, chord, n, spacing, units, outdir, p_fun, wall, nInner, pscale)
%AIRFOIL_DISCRETIZE  Fit ONE continuous spline through an airfoil .dat file
% (Selig or Lednicer format), split it into N arcs of controlled length, and
% export per-segment geometry, outward normals and lengths.
%
% Two features are created in the same run and written to the same DXF:
%   1. the outer outline, as N SPLINE entities (layer SEGMENTS)
%   2. the inner wall, offset WALL (default 0.3) inward from the outline, as
%      ONE continuous closed polyline (layer INNER)
% so the imported profile is a shell of constant wall thickness. Where the
% airfoil is thinner than twice the wall (near the trailing edge) the upper
% and lower offsets meet; the inner polyline ends in a corner there and the
% section aft of it is solid.
%
%   seg = airfoil_discretize('e423.dat', 25.5, 40, 'cosine', 'in', 'out');
%
%   With pressure data it also graphs the pressure on every segment. The
%   easiest way is to give the XFOIL Cp file (CPWR output) as the 7th input
%   and the Cp-to-pressure multiplier n*K*q as the 10th:
%   airfoil_discretize('e423.dat', 25.5, 40, 'cosine', 'in', 'out', ...
%                      '../e423_a12.cp', 0.3, 200, 16993.2);
%   (16993.2 Pa = K * q = 24 * 0.5 * 1.225 * 34^2; leave it out to get the
%   pressures as plain Cp.)
%
%   The 7th input can instead be your own pressure model. Either way the
%   average pressure on each segment and the net force are returned:
%   [seg, pbar, F] = airfoil_discretize('e423.dat', 25.5, 40, 'cosine', 'in', ...
%                                       'out', @my_pressure);
%   where  p = my_pressure(s, xy, isUpper)  receives column vector s (spline
%   parameter), Nx2 points xy (CAD units) and logical isUpper (true on the
%   upper surface), and returns an Nx1 pressure vector.
%
%   datfile : airfoil file (does not have to be normalized or start at x = 0)
%   chord   : chord length in CAD units, default 25.5. The leading edge
%             (minimum-x point of the spline) is moved to the origin (0, 0)
%             and the outline is scaled so the chord, measured along x from
%             the leading edge to the trailing edge, equals this value.
%   n       : number of segments
%   spacing : 'cosine' (clustered at LE and TE on each surface) or 'uniform'
%   units   : 'in' or 'mm' (DXF header setting only)
%   outdir  : output folder
%   p_fun   : (optional) pressure input: an XFOIL Cp file name ('x Cp' or
%             'x y Cp' rows, same airfoil), or a function handle as above,
%             or [] for none
%   wall    : wall thickness in CAD units, default 0.3 (0 = no inner wall)
%   nInner  : number of vertices in the inner wall polyline, default 200
%   pscale  : used with a Cp file only: pressure = pscale * Cp. Give n*K*q
%             in Pa (load factor * rib spacing/thickness * 0.5*rho*V^2).
%             Default 1, which leaves the pressures as Cp.
%
%   Fourth output inner is the K-by-2 list of inner wall vertices (closed:
%   the last vertex joins back to the first).
%
% Outputs (in outdir):
%   segments.csv        per-segment table (length, midpoint, normal, tangent)
%   inner_polyline.csv  x, y vertices of the inner wall
%   outline.dxf         one SPLINE entity per segment (+ TE closing line if
%                       blunt) + one closed LWPOLYLINE for the inner wall
%   dm_curves.txt       DesignModeler "3D Curve" coordinate file: group,id,x,y,z
%                       (outer segments, then the inner wall as the last group)
%   outline.png       plot for checking breakpoints and normals
%   summary.png       graph of the relative geometry (x/c, y/c: outline,
%                     break points, inner wall, thickness and camber) and,
%                     when pressure is given, the pressure on the segments
%                     (arrows on the airfoil, and average vs local values
%                     along the chord)
%
% Base MATLAB only (no toolboxes). Also runs in GNU Octave.

if nargin < 2 || isempty(chord),   chord = 25.5;       end
if nargin < 3 || isempty(n),       n = 40;             end
if nargin < 4 || isempty(spacing), spacing = 'cosine'; end
if nargin < 5 || isempty(units),   units = 'in';       end
if nargin < 6 || isempty(outdir),  outdir = 'out';     end
if nargin < 7, p_fun = []; end
if nargin < 8 || isempty(wall),    wall = 0.3;         end
if nargin < 9 || isempty(nInner),  nInner = 200;       end
if nargin < 10 || isempty(pscale), pscale = 1;         end
if ~exist(outdir, 'dir'), mkdir(outdir); end
pbar = [];  F = [];

if ~(chord > 0), error('chord must be a positive length.'); end
if n < 2, error('n must be at least 2 (one segment per surface).'); end
if wall < 0, error('wall must be positive (or 0 for no inner wall).'); end
if wall > 0 && nInner < 8, error('nInner must be at least 8.'); end

% Selig order; leading edge -> origin, chord (LE to TE along x) -> chord
[pts, c0] = place_on_chord(read_dat(datfile), chord);
[pp, smax] = build_spline(pts);
[sEdges, Ltot] = breakpoints(pp, smax, n, spacing);
seg = segment_data(pp, sEdges, 25);
sLE = sEdges(floor(n/2) + 1);                      % LE breakpoint

punit = '';  ploc = {};
if ischar(p_fun)                                   % XFOIL Cp file
    if ~(pscale > 0), error('pscale must be a positive number.'); end
    sNodes = [0; cumsum(sqrt(sum(diff(pts, 1, 1).^2, 2)))];
    cpNodes = read_cp(p_fun, pts, chord);
    p_fun = @(s, xy, isUpper) pscale * interp1(sNodes, cpNodes, s);
    if pscale == 1, punit = 'Cp'; else, punit = 'Pa'; end
    fprintf('pressure from %d Cp values, p = %g * Cp\n', numel(cpNodes), pscale);
elseif ~isempty(p_fun)
    punit = 'pressure';
end
if ~isempty(p_fun)
    pbar = average_pressure(pp, seg, p_fun, sLE);
    F    = resultant_force(seg, pbar);
    fprintf('net force per unit span: Fx = %.5g, Fy = %.5g  (pressure x %s)\n', F(1), F(2), units);
    if ~strcmp(punit, 'pressure')
        fprintf(['force coefficients from the segment averages: Cn = %+.4f, Ca = %+.4f ' ...
                 '(compare with XFOIL / cp_visualizer)\n'], F(2)/(pscale*chord), F(1)/(pscale*chord));
    end
    ploc = cell(1, n);                             % local (not averaged) pressure
    for k = 1:n
        sk = linspace(seg(k).s0, seg(k).s1, size(seg(k).xy, 1)).';
        pk = p_fun(sk, seg(k).xy, sk <= sLE);
        ploc{k} = pk(:);
    end
end

inner = [];
if wall > 0
    [inner, info] = inner_wall(pp, smax, pts, wall, nInner);
    write_inner_csv(inner, fullfile(outdir, 'inner_polyline.csv'));
end

write_csv(seg, fullfile(outdir, 'segments.csv'), pbar);
write_dm_curves(seg, fullfile(outdir, 'dm_curves.txt'), inner);
write_dxf(seg, pts, fullfile(outdir, 'outline.dxf'), units, inner);
try
    make_plot(seg, pts, fullfile(outdir, 'outline.png'), inner);
catch err
    warning('Plot skipped: %s', err.message);
end
try
    make_summary(seg, inner, chord, units, wall, pbar, ploc, punit, fullfile(outdir, 'summary.png'));
catch err
    warning('Summary graph skipped: %s', err.message);
end

le = seg(floor(n/2) + 1).xy(1, :);
te = 0.5 * (pts(1,:) + pts(end,:));
fprintf('input chord: %.5f  -> scaled by %.5f to chord %g %s\n', c0, chord / c0, chord, units);
fprintf('LE at (%.6f, %.6f)  | TE at (%.6f, %.6f)\n', le(1), le(2), te(1), te(2));
L = [seg.length];
fprintf('points read: %d  | outline arc length: %.4f\n', size(pts,1), Ltot);
fprintf('segments: %d  | length min/max: %.4f / %.4f\n', n, min(L), max(L));
fprintf('TE gap (first-last point): %.5f\n', norm(pts(1,:) - pts(end,:)));
if ~isempty(inner)
    fprintf(['inner wall: %g %s inside the outline, one closed polyline of %d ' ...
             'vertices, length %.4f\n'], wall, units, size(inner,1), info.length);
    fprintf(['inner wall spans x = %.4f to %.4f (solid aft of that) | thinnest ' ...
             'wall at a vertex: %.4f\n'], info.x_le, info.x_te, info.min_gap);
    for k = 1:size(info.mid_corners, 1)
        fprintf(['[NOTE] the outline curves tighter than the %g wall near (%.3f, %.3f); ' ...
                 'the inner wall has a sharp corner there instead of a rounded offset ' ...
                 '(the wall is slightly thicker than nominal at that spot).\n'], ...
                 wall, info.mid_corners(k,1), info.mid_corners(k,2));
    end
end
end


%% ------------------------------------------------------------- reading
function pts = read_dat(path)
% Returns Nx2 array ordered TE(upper) -> LE -> TE(lower) (Selig order).
fid = fopen(path, 'r');
if fid < 0, error('Cannot open %s', path); end
rows = zeros(0, 2);
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    ln = strrep(ln, ',', ' ');
    v = sscanf(ln, '%f');          % non-numeric header lines give []
    if numel(v) == 2 && isempty(regexp(strtrim(ln), '[A-Za-z]', 'once'))
        rows(end+1, :) = v(:).'; %#ok<AGROW>
    end
end
fclose(fid);
pts = rows;

% Lednicer: first numeric row is the point counts (e.g. "31. 31."), then
% upper LE->TE, then lower LE->TE.
if pts(1,1) > 1.5 && pts(1,2) > 1.5
    nu = round(pts(1,1));  nl = round(pts(1,2));
    upper = pts(2:1+nu, :);              % LE -> TE
    lower = pts(2+nu:1+nu+nl, :);        % LE -> TE
    pts = [flipud(upper); lower(2:end, :)];   % TE(up) -> LE -> TE(low)
end
end


function cp = read_cp(path, pts, chord)
% Read an XFOIL Cp file ('x Cp' or 'x y Cp' rows; header and comment lines
% are skipped) and return one Cp per airfoil point, in the airfoil's order.
% pts are the airfoil points already placed on the chord (LE at origin).
fid = fopen(path, 'r');
if fid < 0, error('Cannot open Cp file %s', path); end
rows = zeros(0, 2);
while true
    ln = fgetl(fid);
    if ~ischar(ln), break; end
    ln = strtrim(strrep(ln, ',', ' '));
    if isempty(ln) || ~isempty(regexp(ln, '^[#!%]', 'once')) || ...
            ~isempty(regexp(ln, '[A-CF-Za-cf-z]', 'once'))
        continue;                              % comment or header line
    end
    v = sscanf(ln, '%f');
    if numel(v) == 2 || numel(v) == 3
        rows(end+1, :) = [v(1), v(end)]; %#ok<AGROW>
    end
end
fclose(fid);
if size(rows, 1) < 10
    error('Cp file %s has only %d usable rows (need x Cp or x y Cp columns).', path, size(rows, 1));
end
xc = rows(:,1);  cpv = rows(:,2);
xa = pts(:,1) / chord;                         % airfoil x/c, LE at 0
xc = (xc - min(xc)) / (max(xc) - min(xc));     % Cp file x/c, LE at 0
tol = 2e-3;
if numel(xc) == numel(xa) && max(abs(xc - xa)) < tol
    cp = cpv;
elseif numel(xc) == numel(xa) && max(abs(flipud(xc) - xa)) < tol
    cp = flipud(cpv);
    warning('Cp file is ordered opposite to the airfoil points; reversed it.');
else
    % different point set: split each at its own leading edge and
    % interpolate in x/c on the upper and lower surfaces separately
    warning(['Cp file has %d points and the airfoil %d, or their x stations differ; ' ...
             'Cp was interpolated in x/c per surface. Check summary.png.'], numel(xc), numel(xa));
    [~, ic] = min(xc);  [~, ia] = min(xa);
    if mean(cpv(1:ic)) > mean(cpv(ic:end))     % make the Cp file upper-first
        xc = flipud(xc);  cpv = flipud(cpv);  [~, ic] = min(xc);
    end
    cp = zeros(size(xa));
    [xu, iu] = unique(xc(1:ic));    cu = cpv(1:ic);    cu = cu(iu);
    [xl, il] = unique(xc(ic:end));  cl = cpv(ic:end);  cl = cl(il);
    cp(1:ia)   = interp1(xu, cu, min(max(xa(1:ia), xu(1)), xu(end)));
    cp(ia:end) = interp1(xl, cl, min(max(xa(ia:end), xl(1)), xl(end)));
end
end


%% ------------------------------------------------------------ geometry
function [pp, smax] = build_spline(pts)
% One continuous cubic spline [x(s); y(s)], s = cumulative chord length.
% (MATLAB's spline() uses not-a-knot end conditions.)
d = sqrt(sum(diff(pts, 1, 1).^2, 2));
s = [0; cumsum(d)];
pp = spline(s, pts.');          % 2-by-N data -> vector-valued piecewise poly
smax = s(end);
end

function dpp = ppder1(pp)
% Derivative of a vector-valued pp (base MATLAB; fnder needs a toolbox).
[b, c, l, k, dd] = unmkpp(pp);
dc = c(:, 1:k-1) .* repmat(k-1:-1:1, size(c, 1), 1);
dpp = mkpp(b, dc, dd);
end

function [s, L] = arclength_table(pp, smax)
s = linspace(0, smax, 20000).';
xy = ppval(pp, s).';
ds = sqrt(sum(diff(xy, 1, 1).^2, 2));
L = [0; cumsum(ds)];
end

function sLE = find_le(pp, sTab)
% Spline parameter of the leading edge = minimum-x point of the spline (not
% just the nearest data point): coarse search on the dense table, then a
% bounded 1-D minimisation between the neighbouring table points.
xy = ppval(pp, sTab).';
[xmin, i] = min(xy(:,1));
lo = sTab(max(i-1, 1));  hi = sTab(min(i+1, numel(sTab)));
xfun = @(s) [1 0] * ppval(pp, s);
sLE = fminbnd(xfun, lo, hi, optimset('TolX', 1e-13 * max(sTab(end), 1)));
if xfun(sLE) > xmin, sLE = sTab(i); end
end

function [pts, c0] = place_on_chord(pts, chord)
% Move the leading edge to the origin and scale uniformly so the chord (x
% distance from the leading edge to the trailing edge) equals chord. The
% trailing edge is the midpoint of the first and last points, so a blunt
% trailing edge is handled too. No rotation is applied.
[pp, smax] = build_spline(pts);
sTab = arclength_table(pp, smax);
le = ppval(pp, find_le(pp, sTab)).';
te = 0.5 * (pts(1,:) + pts(end,:));
c0 = te(1) - le(1);
if ~(c0 > 0)
    error(['Trailing edge is not downstream (+x) of the leading edge; ' ...
           'the file must run TE -> upper -> LE -> lower -> TE.']);
end
pts = (pts - repmat(le, size(pts,1), 1)) * (chord / c0);
end

function [sEdges, Ltot] = breakpoints(pp, smax, n, spacing)
% n+1 breakpoints (as spline parameter values). A breakpoint is always put
% exactly at the leading edge (min x), so no segment straddles the LE and
% upper/lower stay separate.
%   uniform : equal true arc length on each surface
%   cosine  : arc-length cosine spacing on EACH surface -> smallest segments
%             at the leading and trailing edges, where dp/ds is steepest.
[sTab, LTab] = arclength_table(pp, smax);
sLE = find_le(pp, sTab);
Lle = interp1(sTab, LTab, sLE);  Ltot = LTab(end);
nUp = floor(n/2);  nLo = n - nUp;

Lup = spaced(nUp, spacing) * Lle;                       % TE -> LE
Llo = Lle + spaced(nLo, spacing) * (Ltot - Lle);        % LE -> TE
target = [Lup; Llo(2:end)];
sEdges = interp1(LTab, sTab, target);
sEdges(nUp + 1) = sLE;              % breakpoint exactly on the leading edge
end

function u = spaced(m, spacing)
u = linspace(0, 1, m + 1).';
if strcmp(spacing, 'cosine')
    u = 0.5 * (1 - cos(pi * u));
end
end

function seg = segment_data(pp, sEdges, npts)
dpp = ppder1(pp);
n = numel(sEdges) - 1;
seg = struct('i', cell(1, n));
for i = 1:n
    s  = linspace(sEdges(i), sEdges(i+1), npts).';
    xy = ppval(pp, s).';
    len = sum(sqrt(sum(diff(xy, 1, 1).^2, 2)));
    sm  = 0.5 * (sEdges(i) + sEdges(i+1));
    mid = ppval(pp, sm).';
    t   = ppval(dpp, sm).';
    t   = t / norm(t);
    % Selig order is counter-clockwise for chord along +x, so the outward
    % normal is the tangent rotated -90 deg: (ty, -tx).
    nv  = [t(2), -t(1)];
    seg(i).i = i - 1;   seg(i).s0 = sEdges(i);  seg(i).s1 = sEdges(i+1);
    seg(i).length = len;  seg(i).xy = xy;
    seg(i).mid = mid;  seg(i).t = t;  seg(i).n = nv;
end
end


%% --------------------------------------------------------- inner wall
function d = dist_to_polyline(Q, A, B)
% Shortest distance from each row of Q to the segments A(i,:) -> B(i,:).
AB = B - A;
den = max(sum(AB.^2, 2), realmin);
d = zeros(size(Q,1), 1);
for k = 1:size(Q,1)
    qa = [Q(k,1) - A(:,1), Q(k,2) - A(:,2)];
    u = min(max(sum(qa .* AB, 2) ./ den, 0), 1);
    r = qa - [u .* AB(:,1), u .* AB(:,2)];
    d(k) = sqrt(min(sum(r.^2, 2)));
end
end

function c = seg_intersection(p1, p2, p3, p4)
% Intersection point of segments p1-p2 and p3-p4, or [] if they do not cross.
c = [];
r = p2 - p1;  q = p4 - p3;
den = r(1)*q(2) - r(2)*q(1);
if abs(den) < realmin, return; end
w = p3 - p1;
a = (w(1)*q(2) - w(2)*q(1)) / den;
b = (w(1)*r(2) - w(2)*r(1)) / den;
if a >= -1e-9 && a <= 1 + 1e-9 && b >= -1e-9 && b <= 1 + 1e-9
    c = p1 + a * r;
end
end

function [V, info] = inner_wall(pp, smax, pts, wall, nInner)
% Inner wall = the outline offset WALL inward, as ONE closed polyline.
%  1. Sample the outer spline densely and step each point WALL along the
%     inward normal (raw offset curve).
%  2. Drop every raw point that is closer than WALL to ANY part of the
%     outline (including the closing line of a blunt trailing edge). This
%     removes the self-crossing loop near the trailing edge, and anywhere
%     else the outline is too thin or too tightly curved for the wall.
%  3. Where a stretch was dropped, join the two remaining ends at the point
%     where the raw offset crosses itself (a sharp corner).
%  4. Resample to nInner vertices, clustered at the inner leading edge and
%     at the trailing-edge corner like the outer segments.
m = 4000;
s = linspace(0, smax, m).';
P = ppval(pp, s).';
T = ppval(ppder1(pp), s).';
T = T ./ repmat(sqrt(sum(T.^2, 2)), 1, 2);
Q = P - wall * [T(:,2), -T(:,1)];             % outward n = (ty, -tx)

A = P(1:end-1, :);  B = P(2:end, :);
if norm(pts(1,:) - pts(end,:)) > 1e-9          % blunt TE closing line
    A = [A; P(end,:)];  B = [B; P(1,:)];
end
idx = find(dist_to_polyline(Q, A, B) >= wall * (1 - 1e-3));
if numel(idx) < 10
    error(['wall = %g is too thick for this airfoil at this chord: no inner ' ...
           'wall is left. Reduce wall.'], wall);
end

loop = zeros(0, 2);  corners = zeros(0, 2);
kf = idx(1);  kl = idx(end);
if kf > 1 && kl < m                            % corner at the trailing edge
    loop = [loop; seg_intersection(Q(kl,:), Q(kl+1,:), Q(kf-1,:), Q(kf,:))];
end
for k = 1:numel(idx) - 1
    a = idx(k);  b = idx(k+1);
    loop(end+1, :) = Q(a,:); %#ok<AGROW>
    if b > a + 1                               % dropped stretch mid-outline
        c = seg_intersection(Q(a,:), Q(a+1,:), Q(b-1,:), Q(b,:));
        loop = [loop; c]; %#ok<AGROW>
        if isempty(c), c = Q(a,:); end
        corners(end+1, :) = c; %#ok<AGROW>
    end
end
loop(end+1, :) = Q(kl,:);

% resample: start vertex (TE end) -> inner LE -> back to start
closed = [loop; loop(1,:)];
L = [0; cumsum(sqrt(sum(diff(closed, 1, 1).^2, 2)))];
keep = [true; diff(L) > 0];                    % interp1 needs distinct L
L = L(keep);  closed = closed(keep, :);
[xle, iLE] = min(closed(:,1));
Lle = L(iLE);
nUp = floor(nInner / 2);
up = spaced(nUp, 'cosine') * Lle;
lo = Lle + spaced(nInner - nUp, 'cosine') * (L(end) - Lle);
target = [up; lo(2:end-1)];
V = [interp1(L, closed(:,1), target), interp1(L, closed(:,2), target)];

info.x_te = loop(1,1);  info.x_le = xle;  info.length = L(end);
info.mid_corners = corners;
info.min_gap = min(dist_to_polyline(V, A, B));
end


%% ----------------------------------------------------------- pressure
function pbar = average_pressure(pp, seg, p_fun, sLE, nint)
% p_bar_i = (1/L_i) * integral p ds over segment i (trapezoid rule).
% p_fun(s, xy, isUpper) -> pressure at those points (see header).
if nargin < 5, nint = 200; end
pbar = zeros(numel(seg), 1);
for i = 1:numel(seg)
    s  = linspace(seg(i).s0, seg(i).s1, nint).';
    xy = ppval(pp, s).';
    ds = sqrt(sum(diff(xy, 1, 1).^2, 2));
    p  = p_fun(s, xy, s <= sLE);
    p  = p(:);
    pbar(i) = sum(0.5 * (p(2:end) + p(1:end-1)) .* ds) / sum(ds);
end
end

function F = resultant_force(seg, pbar)
% Force on the body per unit span from pressure acting INTO the surface,
% n outward:  F_i = -p_i * INT(n ds).  For an arc, INT(n ds) = (dy, -dx)
% with (dx,dy) = end - start point, so this is exact for constant pressure
% on a curved segment (what Mechanical sees with a constant "Normal To"
% pressure). Compare the total to your expected lift and drag.
F = [0 0];
for i = 1:numel(seg)
    d = seg(i).xy(end, :) - seg(i).xy(1, :);
    F = F - pbar(i) * [d(2), -d(1)];
end
end


%% ------------------------------------------------------------- output
function write_csv(seg, path, pbar)
fid = fopen(path, 'w');
hdr = 'i,s0,s1,length,x_mid,y_mid,nx,ny,tx,ty';
if ~isempty(pbar), hdr = [hdr ',p_avg']; end
fprintf(fid, '%s\n', hdr);
for k = 1:numel(seg)
    g = seg(k);
    fprintf(fid, '%d,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f,%.6f', ...
        g.i, g.s0, g.s1, g.length, g.mid(1), g.mid(2), g.n(1), g.n(2), g.t(1), g.t(2));
    if ~isempty(pbar), fprintf(fid, ',%.8g', pbar(k)); end
    fprintf(fid, '\n');
end
fclose(fid);
end

function write_inner_csv(inner, path)
fid = fopen(path, 'w');
fprintf(fid, 'i,x,y\n');
for k = 1:size(inner, 1)
    fprintf(fid, '%d,%.6f,%.6f\n', k - 1, inner(k,1), inner(k,2));
end
fclose(fid);
end

function write_dm_curves(seg, path, inner)
% DesignModeler "3D Curve" coordinate file: Curve Group, Point ID, X, Y, Z.
% One group per segment so each arc is its own edge. The inner wall is the
% last group; its final line has point ID 0, which closes the curve.
fid = fopen(path, 'w');
for k = 1:numel(seg)
    for j = 1:size(seg(k).xy, 1)
        fprintf(fid, '%d\t%d\t%.8f\t%.8f\t%.8f\n', k, j, seg(k).xy(j,1), seg(k).xy(j,2), 0);
    end
end
if ~isempty(inner)
    g = numel(seg) + 1;
    for j = 1:size(inner, 1)
        fprintf(fid, '%d\t%d\t%.8f\t%.8f\t%.8f\n', g, j, inner(j,1), inner(j,2), 0);
    end
    fprintf(fid, '%d\t%d\n', g, 0);
end
fclose(fid);
end

function write_dxf(seg, pts, path, units, inner)
% Minimal ASCII DXF (R2000) with one fit-point SPLINE per segment.
% Hand-written: it has been checked by loading it in a DXF reader, but not
% in Ansys. If Ansys rejects it, use dm_curves.txt instead.
code = 1; if strcmp(units, 'mm'), code = 4; end
fid = fopen(path, 'w');
w = @(c, v) fprintf(fid, '%d\n%s\n', c, num2str(v, '%.10g'));
ws = @(c, v) fprintf(fid, '%d\n%s\n', c, v);

ws(0, 'SECTION'); ws(2, 'HEADER');
ws(9, '$ACADVER'); ws(1, 'AC1015');
ws(9, '$INSUNITS'); w(70, code);
ws(0, 'ENDSEC');
ws(0, 'SECTION'); ws(2, 'ENTITIES');
h = 256;
for k = 1:numel(seg)
    xy = seg(k).xy;  nfit = size(xy, 1);
    ws(0, 'SPLINE'); ws(5, dec2hex(h)); h = h + 1;
    ws(100, 'AcDbEntity'); ws(8, 'SEGMENTS');
    ws(100, 'AcDbSpline');
    w(210, 0); w(220, 0); w(230, 1);
    w(70, 1064);              % planar + fit-point spline
    w(71, 3);                 % degree
    w(72, 0); w(73, 0);       % knots, control points (none: fit-point form)
    w(74, nfit);
    w(42, 1e-7); w(43, 1e-7); w(44, 1e-10);
    for j = 1:nfit
        w(11, xy(j,1)); w(21, xy(j,2)); w(31, 0);
    end
end
if norm(pts(1,:) - pts(end,:)) > 1e-9      % blunt TE: close with a line
    ws(0, 'LINE'); ws(5, dec2hex(h)); ws(100, 'AcDbEntity'); ws(8, 'SEGMENTS');
    ws(100, 'AcDbLine');
    w(10, pts(end,1)); w(20, pts(end,2)); w(30, 0);
    w(11, pts(1,1));   w(21, pts(1,2));   w(31, 0);
    h = h + 1;
end
if ~isempty(inner)                         % inner wall: ONE closed polyline
    ws(0, 'LWPOLYLINE'); ws(5, dec2hex(h)); ws(100, 'AcDbEntity'); ws(8, 'INNER');
    ws(100, 'AcDbPolyline');
    w(90, size(inner, 1));    % vertex count
    w(70, 1);                 % closed
    for j = 1:size(inner, 1)
        w(10, inner(j,1)); w(20, inner(j,2));
    end
end
ws(0, 'ENDSEC'); ws(0, 'EOF');
fclose(fid);
end

function make_plot(seg, pts, path, inner)
f = figure('Visible', 'off', 'Position', [100 100 900 320]);
hold on;
plot(pts(:,1), pts(:,2), 'k.', 'MarkerSize', 4);
c = max(pts(:,1)) - min(pts(:,1));
for k = 1:numel(seg)
    plot(seg(k).xy(:,1), seg(k).xy(:,2), '-', 'LineWidth', 1.2);
    plot(seg(k).xy(1,1), seg(k).xy(1,2), 'r+', 'MarkerSize', 6);
    quiver(seg(k).mid(1), seg(k).mid(2), 0.04*c*seg(k).n(1), 0.04*c*seg(k).n(2), 0, 'g');
end
ttl = sprintf('%d segments (red = breakpoints, green = outward normals)', numel(seg));
if ~isempty(inner)
    plot([inner(:,1); inner(1,1)], [inner(:,2); inner(1,2)], 'k-', 'LineWidth', 0.9);
    ttl = sprintf('%s + inner wall polyline (%d vertices)', ttl, size(inner, 1));
end
axis equal;
title(ttl);
print(f, path, '-dpng', '-r150');
close(f);
end

function make_summary(seg, inner, chord, units, wall, pbar, ploc, punit, path)
% One figure: the geometry in relative coordinates (everything divided by
% the chord, leading edge at the origin) and, when pressure is given, the
% pressure on the segments.
n = numel(seg);  nUp = floor(n/2);
hasP = ~isempty(pbar);
C_SUC = [0.75 0.22 0.17];  C_PRE = [0.12 0.37 0.66];  C_BODY = [0.17 0.17 0.17];

% dense outline and segment end points, relative to the chord
XY = zeros(0, 2);
for k = 1:n, XY = [XY; seg(k).xy / chord]; end %#ok<AGROW>
brk = zeros(n + 1, 2);
for k = 1:n, brk(k, :) = seg(k).xy(1, :) / chord; end
brk(n + 1, :) = seg(n).xy(end, :) / chord;

% thickness and camber from the upper and lower surfaces at common x/c
up = zeros(0, 2);  lo = zeros(0, 2);
for k = 1:nUp,   up = [up; seg(k).xy / chord]; end %#ok<AGROW>
for k = nUp+1:n, lo = [lo; seg(k).xy / chord]; end %#ok<AGROW>
[xu, iu] = unique(up(:,1));  yu = up(iu, 2);
[xl, il] = unique(lo(:,1));  yl = lo(il, 2);
xq = linspace(max(xu(1), xl(1)), min(xu(end), xl(end)), 400).';
yuq = interp1(xu, yu, xq);  ylq = interp1(xl, yl, xq);
thick = yuq - ylq;  camber = 0.5 * (yuq + ylq);
[tmax, it] = max(thick);  [cmax, ic] = max(camber);

if hasP, rows = 3; hgt = 900; else, rows = 1; hgt = 360; end
f = figure('Visible', 'off', 'Position', [100 100 1000 hgt]);

% ---- 1. relative geometry
subplot(rows, 1, 1);  hold on;
plot(XY(:,1), XY(:,2), '-', 'Color', C_BODY, 'LineWidth', 1.3);
plot(xq, camber, '--', 'Color', [0.45 0.45 0.45], 'LineWidth', 0.8);
if ~isempty(inner)
    plot([inner(:,1); inner(1,1)] / chord, [inner(:,2); inner(1,2)] / chord, '-', ...
         'Color', [0.20 0.55 0.30], 'LineWidth', 1.0);
end
plot(brk(:,1), brk(:,2), '.', 'Color', C_SUC, 'MarkerSize', 9);
plot([xq(it) xq(it)], [ylq(it) yuq(it)], '-', 'Color', C_PRE, 'LineWidth', 1.2);
plot(0, 0, 'ko', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
plot(brk(1,1), brk(1,2), 'ks', 'MarkerSize', 6, 'MarkerFaceColor', 'k');
text(0.01, -0.035, 'LE (0, 0)', 'FontSize', 8);
text(brk(1,1) - 0.01, brk(1,2) - 0.035, 'TE', 'FontSize', 8, 'HorizontalAlignment', 'right');
text(xq(it) + 0.01, ylq(it) - 0.025, sprintf('max thickness %.1f%% c at x/c = %.2f', ...
     100 * tmax, xq(it)), 'FontSize', 8, 'Color', C_PRE);
grid on;
set(gca, 'XLim', [-0.05 1.05], 'YLim', [min(XY(:,2)) - 0.06, max(XY(:,2)) + 0.03], ...
         'DataAspectRatio', [1 1 1]);
xlabel('x / c');  ylabel('y / c');
ttl = sprintf(['Relative geometry: chord c = %g %s, %d segments (dots = break points), ' ...
               'max camber %.1f%% c at x/c = %.2f'], chord, units, n, 100 * cmax, xq(ic));
if ~isempty(inner)
    ttl = sprintf('%s, inner wall %.4f c (%g %s)', ttl, wall / chord, wall, units);
end
title(ttl, 'FontSize', 9);

if hasP
    pmax = max(abs(pbar));  if pmax == 0, pmax = 1; end
    if strcmp(punit, 'Cp'), plabel = 'Cp'; else, plabel = sprintf('pressure [%s]', punit); end

    % ---- 2. average pressure on each segment, drawn on the airfoil
    % p < 0 (suction) pulls OUTWARD: red arrow pointing away from the surface.
    % p > 0 pushes INTO the surface: blue arrow pointing at the surface.
    subplot(rows, 1, 2);  hold on;
    plot(XY(:,1), XY(:,2), '-', 'Color', C_BODY, 'LineWidth', 1.3);
    scale = 0.25 / pmax;                           % longest arrow = 0.25 chord
    tips = XY;
    for k = 1:n
        m = seg(k).mid / chord;  len = abs(pbar(k)) * scale;  nv = seg(k).n;
        far = m + len * nv;                        % arrow end away from the surface
        tips(end+1, :) = far; %#ok<AGROW>
        if pbar(k) < 0
            draw_arrow(m, far, C_SUC);             % tail on the surface, head outward
        else
            draw_arrow(far, m, C_PRE);             % head on the surface
        end
    end
    plot(brk(:,1), brk(:,2), '.', 'Color', C_BODY, 'MarkerSize', 6);
    grid on;
    set(gca, 'XLim', [min(tips(:,1)) - 0.05, max(tips(:,1)) + 0.05], ...
             'YLim', [min(tips(:,2)) - 0.03, max(tips(:,2)) + 0.03], ...
             'DataAspectRatio', [1 1 1]);
    xlabel('x / c');  ylabel('y / c');
    title(sprintf(['Average %s on each segment (red = suction, pulls outward; blue = pushes ' ...
                   'into the surface); longest arrow = %.4g'], plabel, pmax), 'FontSize', 9);

    % ---- 3. average vs local pressure along the chord
    subplot(rows, 1, 3);  hold on;
    hl = zeros(1, 4);
    for k = 1:n
        if k <= nUp, col = C_SUC; else, col = C_PRE; end
        xk = seg(k).xy(:,1) / chord;
        h1 = plot(xk, ploc{k}, '-', 'Color', col, 'LineWidth', 0.8);
        h2 = plot([xk(1) xk(end)], [pbar(k) pbar(k)], '-', 'Color', col, 'LineWidth', 2.4);
        if k == 1,       hl(1) = h1;  hl(2) = h2; end
        if k == nUp + 1, hl(3) = h1;  hl(4) = h2; end
    end
    plot([0 1], [0 0], 'k:');
    set(gca, 'YDir', 'reverse');                   % suction upward, as for Cp plots
    grid on;  xlim([-0.05 1.05]);
    xlabel('x / c');  ylabel(plabel);
    legend(hl, {'upper, local', 'upper, segment average', 'lower, local', ...
                'lower, segment average'}, 'Location', 'northeast');
    title('Local (thin) and segment-average (thick) pressure; suction is up', 'FontSize', 9);
end

print(f, path, '-dpng', '-r150');
close(f);
end

function draw_arrow(a, b, col)
% Arrow from a to b drawn with plain lines, so it looks the same in MATLAB
% and Octave. The head is 30% of the arrow length, at most 0.02.
d = b - a;  len = norm(d);
if len <= 0, return; end
u = d / len;  w = [-u(2), u(1)];
h = min(0.3 * len, 0.02);
p1 = b - h * u + 0.4 * h * w;  p2 = b - h * u - 0.4 * h * w;
plot([a(1) b(1)], [a(2) b(2)], '-', 'Color', col, 'LineWidth', 1.0);
plot([p1(1) b(1) p2(1)], [p1(2) b(2) p2(2)], '-', 'Color', col, 'LineWidth', 1.0);
end
