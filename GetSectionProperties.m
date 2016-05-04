function [WingSection] = GetSectionProperties(nacaDigits, chordLength)

%% Setup Primary Structure Contour
numberOfPoints = 200;
airfoil = NacaAirfoil_t(nacaDigits, chordLength, numberOfPoints);

% Truncate reat 20%
tmp = find(airfoil.Geometry.X <= chordLength*0.80);
fprintf('Working points: %d through %d\n', tmp(1), tmp(end));
airfoil.Geometry.Y = airfoil.Geometry.Y(airfoil.Geometry.X <= chordLength*0.80);
airfoil.Geometry.X = airfoil.Geometry.X(airfoil.Geometry.X <= chordLength*0.80);
airfoil.Geometry.X(end+1) = airfoil.Geometry.X(1);
airfoil.Geometry.Y(end+1) = airfoil.Geometry.Y(1);

%% Sking sections, spars, and stringers

% Skin Sections
Sections = struct('p1',[],'p2',[], 'L', [], 'A', [], 'x', [], 'y', []);
t_skin = 0.001016;
for ii = 1:length(airfoil.Geometry.X)-1
    Sections(ii).p1 = [airfoil.Geometry.X(ii) airfoil.Geometry.Y(ii) 0];
    Sections(ii).p2 = [airfoil.Geometry.X(ii+1) airfoil.Geometry.Y(ii+1) 0];
    Sections(ii).L = norm(Sections(ii).p2-Sections(ii).p1);
    Sections(ii).A = t_skin*Sections(ii).L;
    midPoint = (Sections(ii).p1+Sections(ii).p2)/2;
    Sections(ii).x = midPoint(1);
    Sections(ii).y = midPoint(2);
end

% Spars
% Only two spars will be added. One spar will be at the end of the primary
% structure. and the other will be placed such that the structure is less
% succeptible to deformations
findNearest = @(array, value, start) min(abs(array(start:end)-value));
SparLocations = [0.8*chordLength 0.25*chordLength];
t_spar = 0.0025;
Spars = struct('p1',[],'p2',[],'L',[],'A',[],'x',[],'y',[]);
idx = 1;
for ii = 1:length(SparLocations)
    [~, idx] = findNearest(airfoil.Geometry.X, SparLocations(ii),idx);    
    idx = find(airfoil.Geometry.X==airfoil.Geometry.X(idx));
    Spars(ii).p1 = [airfoil.Geometry.X(idx(1)) airfoil.Geometry.Y(idx(1))];
    Spars(ii).p2 = [airfoil.Geometry.X(idx(2)) airfoil.Geometry.Y(idx(2))];
    
    Spars(ii).L = norm(Spars(ii).p2-Spars(ii).p1);
    Spars(ii).A = t_spar*Spars(ii).L;
    
    midPoint = (Spars(ii).p1+Spars(ii).p2)/2;
    Spars(ii).x = midPoint(1); % Must be the same X as p1 and p2
    Spars(ii).y = midPoint(2);
end

% Spar Caps
SparCaps = struct('A',[],'x',[],'y',[]);
SparCapArea = 3e-5; % Spar cap area in m^2
for ii = 1:length(Spars)
    if ii ~= 1
        tmpA = SparCapArea*2;
    else
        tmpA = SparCapArea;        
    end
    % The spar ii will have caps 2*ii-1 and 2*ii e.g Spar#2 has 
    % caps 3 and 4
    capIdx = 2*ii-1; 
    SparCaps(capIdx).x = Spars(ii).p1(1);
    SparCaps(capIdx).y = Spars(ii).p1(2);
    SparCaps(capIdx).A = tmpA;
    SparCaps(capIdx+1).x = Spars(ii).p2(1);
    SparCaps(capIdx+1).y = Spars(ii).p2(2);
    SparCaps(capIdx+1).A = tmpA;
end

% Stringers
% Stringers will be added along the
Stringers = struct('A',[],'x',[],'y',[]);
StringerArea = 5e-6; % Stringer crossectional area in m^2
StringerLocations = airfoil.Geometry.X(20:20:end);
StringerLocations = StringerLocations(1:end-1);
idx = 1;
for ii = 1:length(StringerLocations)
    [~, idxf] = findNearest(airfoil.Geometry.X, StringerLocations(ii),idx);
    idx = idx+idxf-1;
    Stringers(ii).x = airfoil.Geometry.X(idx);
    Stringers(ii).y = airfoil.Geometry.Y(idx);
    Stringers(ii).A = StringerArea;   
end

% Find the centroid of the cell
WingSection.A_sum = sum([Sections(:).A]) + sum([Stringers(:).A]) + sum([Spars(:).A]) + sum([SparCaps(:).A]);
WingSection.xbar = (sum([Sections(:).x].*[Sections(:).A]) + sum([Stringers(:).x].*[Stringers(:).A]) + sum([Spars(:).x].*[Spars(:).A]) + sum([SparCaps(:).x].*[SparCaps(:).A]))/WingSection.A_sum;
WingSection.ybar = (sum([Sections(:).y].*[Sections(:).A]) + sum([Stringers(:).y].*[Stringers(:).A]) + sum([Spars(:).y].*[Spars(:).A]) + sum([SparCaps(:).y].*[SparCaps(:).A]))/WingSection.A_sum;

% figure;
% hold all;
% p1 = plot(airfoil.Geometry.X, airfoil.Geometry.Y, '-o');
% p2 = plot([Stringers(:).x], [Stringers(:).y], 'ro','MarkerSize',5);
% for ii = 1:length(Spars)
%     Spar = Spars(ii);
%     p_spar(ii) = plot([Spar.p1(1) Spar.p2(1)], [Spar.p1(2) Spar.p2(2)], 'g','LineWidth',3);
% end
% p3 = plot([SparCaps(:).x], [SparCaps(:).y], 'ks','MarkerSize',10);
% p4 = plot(xbar,ybar,'rx');
% %p4 = plot([Sections(:).x], [Sections(:).y], 'r.','MarkerSize',5);
% hleg = legend([p1, p2, p_spar(1), p3, p4],'Airfoil Points', 'Stringers', 'Spars', 'Spar Caps', 'Centroid');
% title(airfoil.Name);
% axis equal;
% grid on;
% set(gca,'FontName','Times New Roman','FontSize',15);

%% Shift centroid to the origin

% shift airfoil
airfoil.Geometry.X = airfoil.Geometry.X - WingSection.xbar;
airfoil.Geometry.Y = airfoil.Geometry.Y - WingSection.ybar;

% shift skin sections
for ii = 1:length(Sections)
    Sections(ii).p1 = Sections(ii).p1 - [WingSection.xbar WingSection.ybar 0];
    Sections(ii).p2 = Sections(ii).p2 - [WingSection.xbar WingSection.ybar 0];
    Sections(ii).x = Sections(ii).x - WingSection.xbar;
    Sections(ii).y = Sections(ii).y - WingSection.ybar;
end

% shift spars
for ii = 1:length(Spars)
    Spars(ii).p1 = Spars(ii).p1 - [xbar ybar];
    Spars(ii).p2 = Spars(ii).p2 - [xbar ybar];
    Spars(ii).x = Spars(ii).x - xbar;
    Spars(ii).y = Spars(ii).y - ybar;    
end

% shift spar caps
for ii = 1:length(SparCaps)
    SparCaps(ii).x = SparCaps(ii).x - xbar;
    SparCaps(ii).y = SparCaps(ii).y - ybar;    
end

% shift stringers
for ii = 1:length(Stringers)
    Stringers(ii).x = Stringers(ii).x - xbar;
    Stringers(ii).y = Stringers(ii).y - ybar;    
end

% Find the centroid of the cell
xbar = (sum([Sections(:).x].*[Sections(:).A]) + sum([Stringers(:).x].*[Stringers(:).A]) + sum([Spars(:).x].*[Spars(:).A]) + sum([SparCaps(:).x].*[SparCaps(:).A]))/WingSection.A_sum;
ybar = (sum([Sections(:).y].*[Sections(:).A]) + sum([Stringers(:).y].*[Stringers(:).A]) + sum([Spars(:).y].*[Spars(:).A]) + sum([SparCaps(:).y].*[SparCaps(:).A]))/WingSection.A_sum;

%% Find area moments of inertia

% shift skin sections
for ii = 1:length(Sections)
    Sections(ii).Ixx = Sections(ii).A*Sections(ii).y^2;
    Sections(ii).Iyy = Sections(ii).A*Sections(ii).x^2;
    Sections(ii).Ixy = Sections(ii).A*Sections(ii).x*Sections(ii).y;
end

% shift spars
for ii = 1:length(Spars)
    Spars(ii).Ixx = t_spar*Spars(ii).L^3 + Spars(ii).A*Spars(ii).y^2;
    Spars(ii).Iyy = Spars(ii).A*Spars(ii).x^2;
    Spars(ii).Ixy = Spars(ii).A*Spars(ii).x*Spars(ii).y;
end

% shift spar caps
for ii = 1:length(SparCaps)
    SparCaps(ii).Ixx = SparCaps(ii).A*SparCaps(ii).y^2;
    SparCaps(ii).Iyy = SparCaps(ii).A*SparCaps(ii).x^2;
    SparCaps(ii).Ixy = SparCaps(ii).A*SparCaps(ii).x*SparCaps(ii).y; 
end

% shift stringers
for ii = 1:length(Stringers)
    Stringers(ii).Ixx = Stringers(ii).A*Stringers(ii).y^2;
    Stringers(ii).Iyy = Stringers(ii).A*Stringers(ii).x^2;
    Stringers(ii).Ixy = Stringers(ii).A*Stringers(ii).x*Stringers(ii).y;
end

WingSection.Ixx = sum([Sections(:).Ixx]) + sum([Stringers(:).Ixx]) + sum([Spars(:).Ixx]) + sum([SparCaps(:).Ixx]);
WingSection.Iyy = sum([Sections(:).Iyy]) + sum([Stringers(:).Iyy]) + sum([Spars(:).Iyy]) + sum([SparCaps(:).Iyy]);
WingSection.Ixy = sum([Sections(:).Ixy]) + sum([Stringers(:).Ixy]) + sum([Spars(:).Ixy]) + sum([SparCaps(:).Ixy]);

figure;
hold all;
p1 = plot(airfoil.Geometry.X, airfoil.Geometry.Y, '-o');
p2 = plot([Stringers(:).x], [Stringers(:).y], 'ro','MarkerSize',5);
for ii = 1:length(Spars)
    Spar = Spars(ii);
    p_spar(ii) = plot([Spar.p1(1) Spar.p2(1)], [Spar.p1(2) Spar.p2(2)], 'g','LineWidth',3);
end
p3 = plot([SparCaps(:).x], [SparCaps(:).y], 'ks','MarkerSize',10);
p4 = plot(xbar,ybar,'rx');
%p4 = plot([Sections(:).x], [Sections(:).y], 'r.','MarkerSize',5);
hleg = legend([p1, p2, p_spar(1), p3, p4],'Airfoil Points', 'Stringers', 'Spars', 'Spar Caps', 'Centroid');
title([airfoil.Name ' w/ centroid @(0,0)']);
axis equal;
grid on;
text(-0.5, 0.4, sprintf('Ixx = %f',Ixx));
text(-0.5, 0.37, sprintf('Iyy = %f',Iyy));
text(-0.5, 0.34, sprintf('Ixy = %f',Ixy));
set(gca,'FontName','Times New Roman','FontSize',15);
