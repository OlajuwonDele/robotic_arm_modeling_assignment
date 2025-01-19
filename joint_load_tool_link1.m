clear all
%% User Input when all joint angles = 0
% Mass of each link from SOLIDWORKS
m_6 = 0.08050 ; %kg
m_5 = 0.20129 ; %kg
m_4 = 0.91398; %kg
m_3 = 0.45763; %kg
m_2 = 1.35269 ; %kg
m_1 = 0; %kg

%Applied force f_e and torque t_e acting on the end effector
f_e = [0; 0; 0]; %N
t_e = [0; 0; 0]; %Nm

% Global origin / centre of rotation of the links from SOLIDWORKS
origin_6 = [0; 0; 0] / 1000; %m
origin_5 = [0; 0; 0] / 1000; %m
origin_4 = [0; 0; 0] / 1000; %m
origin_3 = [0; 0; 0] / 1000; %m
origin_2 = [0; 0; 0] / 1000; %m
origin_1 = [0; 0; 0] / 1000; %m

% Global point E (Centre of end effector), 
point_coord = [364.65; -1.57; 461.95] / 1000; %m

% Global mass centre position vector from SOLIDWORKS
mass_cent_6 = [346.05;-1.57; 462.97] / 1000; %m
mass_cent_5 = [291.34; 4.67; 465.32] / 1000; %m
mass_cent_4 = [167.03; 1.60; 470.68] / 1000; %m
mass_cent_3 = [80.27; 12.92; 458.69] / 1000; %m
mass_cent_2 = [68.56; 51.62; 304.90] / 1000; %m
mass_cent_1 = [0; 0; 0] / 1000; %m

% Orientation of each joint
theta_1 = 0; %degrees
theta_2 = 0; %degrees
theta_3 = 0; %degrees
theta_4 = 0; %degrees
theta_5 = 0; %degrees
theta_6 = 0; %degrees

%% Visualization 
%To interpret the configuration of robot’s posture using DH parameters run the following
%section

%a, alpha, d, theta
dhparams = [64.200	-90	169.770	theta_1;
305.000	0.000	0.000	theta_2 - 90;
0.000	90	0.000	theta_3 + 180;
0.000	-90	222.630	theta_4;
0.000	90	0.000	theta_5;
0.000	0.000	41.000	theta_6];

figure(1)
hold on
grid on

t = eye(4);
x = [];
y = [];
z = [];

for i = 1:6

   a = dhparams(i, 1);
   alpha = dhparams(i, 2);
   d = dhparams(i, 3);
   theta = dhparams(i, 4);
   
   current_t = DHConvention(a, alpha, d, theta);
   t = t * current_t;  
   
   position = t(1:3, 4);

   x = [x; position(1)];
   y = [y; position(2)];
   z = [z; position(3)];
   
   scatter3(position(1), position(2), position(3), 'filled');
   plotframe(current_t(1:3, 1:3), position, [30 30 30])

   text(position(1), position(2), position(3), ['Joint ', num2str(i)], 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right', 'Color','red','FontSize', 12);
end


plot3(x, y, z, 'k-', 'LineWidth', 3);
xlim([-500 500])
ylim([-500 500])
zlim([-500 500])
view(3)
xlabel('X');
ylabel('Y');
zlabel('Z');
title('Robot Arm');

hold off
  
%% Joint load and torque calculation
% No need to alter
g = 9.81; % m/s-2 %kg

m_s = 0; %kg

% Joint angle position vector 
% Measure the distance from origin a to origin b from SOLIDWORKS
r_6_E = point_coord - origin_6;
r_5_6 = origin_6 - origin_5;
r_4_5 = origin_5 - origin_4;
r_3_4 = origin_4 - origin_3;
r_2_3 = origin_3 - origin_2;
r_1_2 = origin_2 - origin_1;

R_1 = eye(3);
R_2 = R_1 * rotz(theta_1);
R_3 = R_2 * roty(theta_2);
R_4 = R_3 * roty(theta_3);
R_5 = R_4 * rotx(theta_4);
R_6 = R_5 * roty(theta_5);
R_e = R_6 * rotx(theta_6);

z = [0; 0; 1];

%Rotating coordinate frames and vectors relative to joint angles
r_Ge_E = R_e * (point_coord - mass_cent_6);
r_G5_6 = R_6 * (origin_6 - mass_cent_5);
r_G4_5 = R_5 * (origin_5 - mass_cent_4);
r_G3_4 = R_4 * (origin_4 - mass_cent_3);
r_G2_3 = R_3 * (origin_3 - mass_cent_2);
r_G1_2 = R_2 * (origin_2 - mass_cent_1);

point_coord = R_e * point_coord;
r_6_E = R_e * r_6_E;
r_5_6 = R_6 * r_5_6;
r_4_5 = R_5 * r_4_5;
r_3_4 = R_4 * r_3_4;
r_2_3 = R_3 * r_2_3;
r_1_2 = R_2 * r_1_2;

mass_cent_1 = R_2 * mass_cent_1;
mass_cent_2 = R_3 * mass_cent_2;
mass_cent_3 = R_4 * mass_cent_3;
mass_cent_4 = R_5 * mass_cent_4;
mass_cent_5 = R_6 * mass_cent_5;
mass_cent_6 = R_e * mass_cent_6;

origin_1 = R_2 * origin_1;
origin_2 = R_3 * origin_2;
origin_3 = R_4 * origin_3;
origin_4 = R_5 * origin_4;
origin_5 = R_6 * origin_5;
origin_6 = R_e * origin_6;


%joint 6
m_e = m_s + m_6;
rho_w6f = inv(R_6) * (f_e - ((m_e * g * z)));
rho_w6t = inv(R_6) * (t_e - (cross((m_e * g * (R_e * r_Ge_E)), z) + cross(R_e * r_6_E, R_6 * rho_w6f)));

SwE = [f_e; t_e];
SwG_E = - m_e * g * [z; cross((R_e * r_Ge_E), z)];
Sw6 = SwE + SwG_E;
rho_w6_ = [rho_w6f.', rho_w6t.'].';

r6 = R_e * r_6_E;
rhat_6 = skew(r6);
W_6 = [R_6 rhat_6 * R_6; zeros(3,3) R_6];

rho_w6 = W_6.' * Sw6;

%joint 5
rho_w5f = inv(R_5) * (f_e - (((m_e + m_5) * g * z))) ;
rho_w5t = inv(R_5) *(t_e - (cross(m_e * g * (R_e * r_Ge_E), z) + cross(m_5 * g * (R_e * r_6_E + R_6 * r_G5_6), z)+ cross((R_e * r_6_E) + (R_6 * r_5_6), R_5 * rho_w5f)));

SwE = [f_e; t_e];
SwG_5 = - m_5 * g * [z; cross((r6 + R_6 * r_G5_6), z)];
Sw5 = SwE + SwG_E + SwG_5;
rho_w5_ = [rho_w5f.', rho_w5t.'].';

r5 = r6 + R_6 * r_5_6;
rhat_5 = skew(r5);
W_5 = [R_5 rhat_5 * R_5; zeros(3,3) R_5];
rho_w5 = W_5.' * Sw5;

%joint 4
rho_w4f = inv(R_4) * (f_e - (((m_e + m_5 + m_4) * g * z)));
rho_w4t = inv(R_4) * (t_e - (cross(m_e * g * (R_e * r_Ge_E), z) + cross(m_5 * g * (R_e * r_6_E + R_6 * r_G5_6), z) + cross(m_4 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_G4_5)), z) + cross((R_e * r_6_E) + (R_6 * r_5_6) + (R_5 * r_4_5), R_4 * rho_w4f)));

SwE = [f_e; t_e];
SwG_4 = - m_4 * g * [z; cross((r5 + R_5 * r_G4_5), z)];
Sw4 = SwE + SwG_E + SwG_5 + SwG_4;
rho_w4_ = [rho_w4f.', rho_w4t.'].';

r4 = r5 + R_5 * r_4_5;
rhat_4 = skew(r4);
W_4 = [R_4 rhat_4 * R_4; zeros(3,3) R_4];
rho_w4 = W_4.' * Sw4;


%joint 3
rho_w3f = inv(R_3) *(f_e - (((m_e + m_5 + m_4 + m_3) * g * z)));
rho_w3t = inv(R_3) * (t_e - (cross(m_e * g * (R_e * r_Ge_E), z) + cross(m_5 * g * (R_e * r_6_E + R_6 * r_G5_6), z) + cross(m_4 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_G4_5)), z) + cross(m_3 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_G3_4)), z) + cross((R_e * r_6_E) + (R_6 * r_5_6) + (R_5 * r_4_5) + (R_4 * r_3_4), R_3 * rho_w3f)));

SwE = [f_e; t_e];
SwG_3 = - m_3 * g * [z; cross((r4 + R_4 * r_G3_4), z)];
Sw3 = SwE + SwG_E + SwG_5 + SwG_4 + SwG_3;
rho_w3_ = [rho_w3f.', rho_w3t.'].';

r3 = r4 + R_4 * r_3_4;
rhat_3 = skew(r3);
W_3 = [R_3 rhat_3 * R_3; zeros(3,3) R_3];
rho_w3 = W_3.' * Sw3;

%joint 2
rho_w2f = inv(R_2) * (f_e - (((m_e + m_5 + m_4 + m_3 + m_2) * g * z)));
rho_w2t = inv(R_2) * (t_e - (cross(m_e * g * (R_e * r_Ge_E), z) + cross(m_5 * g * (R_e * r_6_E + R_6 * r_G5_6), z) + cross(m_4 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_G4_5)), z) + cross(m_3 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_G3_4)), z) + cross(m_2 * g * ((R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_3_4 + R_3 * r_G2_3)), z) + cross((R_e * r_6_E) + (R_6 * r_5_6) + (R_5 * r_4_5) + (R_4 * r_3_4) + (R_3 * r_2_3), R_2 * rho_w2f)));

SwE = [f_e; t_e];
SwG_2 = - m_2 * g * [z; cross((r3 + R_3 * r_G2_3), z)];
Sw2 = SwE + SwG_E + SwG_5 + SwG_4 + SwG_3 + SwG_2;
rho_w2 = [rho_w2f.', rho_w2t.'].';

r2 = r3 + R_3 * r_2_3;
rhat_2 = skew(r2);
W_2 = [R_2 rhat_2 * R_2; zeros(3,3) R_2];
rho_w2_ = W_2.' * Sw2;

%joint 1 
rho_w1f = inv(R_1) *(f_e - (((m_e + m_5 + m_4 + m_3 + m_2 + m_1) * g * z)));
rho_w1t = inv(R_1) *(t_e - (cross(m_e * g * (R_e * r_Ge_E), z) + ...
      cross(m_5 * g * (R_e * r_6_E + R_6 * r_G5_6), z) + ...
      cross(m_4 * g * (R_e * r_6_E + R_6 * r_5_6 + R_5 * r_G4_5), z) + ...
      cross(m_3 * g * (R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_G3_4), z) + ...
      cross(m_2 * g * (R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_3_4 + R_3 * r_G2_3), z) + ...
      cross(m_1 * g * (R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_3_4 + R_3 * r_2_3 + R_2 * r_G1_2), z) + ...
      cross(R_e * r_6_E + R_6 * r_5_6 + R_5 * r_4_5 + R_4 * r_3_4 + R_3 * r_2_3 + R_2 * r_1_2, R_1 * rho_w1f) ));

SwE = [f_e; t_e];
SwG_1 = - m_1 * g * [z; cross((r2 + R_2 * r_G1_2), z)];
Sw1 = SwE + SwG_E + SwG_5 + SwG_4 + SwG_3 + SwG_2 + SwG_1;
rho_w1_ = [rho_w1f.', rho_w1t.'].';

r1 = r2 + R_2 * r_1_2;
rhat_1 = skew(r1);
W_1 = [R_1 rhat_1; zeros(3,3) R_1];
rho_w1 = W_1.' * Sw1;

fprintf('The force [N] on Link 1 is: \n');
disp(rho_w1(1:3))

fprintf('The torque [Nm] on Link 1 is: \n');
disp(rho_w1(4:6))

fprintf('The force [N] on Link 2 is: \n');
disp(rho_w2(1:3))

fprintf('The torque [Nm] on Link 2 is: \n');
disp(rho_w2(4:6))

fprintf('The force [N] on Link 3 is: \n');
disp(rho_w3(1:3))

fprintf('The torque [Nm] on Link 3 is: \n');
disp(rho_w3(4:6))

fprintf('The force [N] on Link 4 is: \n');
disp(rho_w4(1:3))

fprintf('The torque [Nm] on Link 4 is: \n');
disp(rho_w4(4:6))

fprintf('The force [N] on Link 5 is: \n');
disp(rho_w5(1:3))

fprintf('The torque [Nm] on Link 5 is: \n');
disp(rho_w5(4:6))

fprintf('The force [N] on Link 6 is: \n');
disp(rho_w6(1:3))

fprintf('The torque [Nm] on Link 6 is: \n');
disp(rho_w6(4:6))

function dh = DHConvention(a,alpha,d,theta)
ca = cosd(alpha);
sa = sind(alpha);

c0 = cosd(theta);
s0 = sind(theta);

dh = [c0 -s0*ca s0*sa a*c0; s0 c0*ca -c0*sa a*s0; 0 sa ca d; 0 0 0 1];
end
%% IGNORE

function S = skew(v)
    if isvec(v,3)
        % SO(3) case
        S = [0   -v(3)  v(2); v(3)  0 -v(1); -v(2) v(1) 0];
    elseif isvec(v,1)
        % SO(2) case
        S = [0 -v; v 0];
    else
        error('SMTB:skew:badarg', 'argument must be a 1- or 3-vector');
    end
end
function h = isvec(v, l)
    if nargin == 1
            l = 3;
    end
    if isa(v, 'symfun')
        h = logical( length(formula(v)) == l);
    else
        d = size(v);
        h = logical( length(d) == 2 && min(d) == 1 && numel(v) == l );
    end
end


function varargout = plotframe( rotationMatrix, translationVector, ...
        basisVectorLengths, Options, QuiverProperties )
%PLOTFRAME Plot a 3-D Cartesian coordinate system.
%     plotframe( )
%     plotframe( rotationMatrix, translationVector )
%     plotframe( rotationMatrix, translationVector, basisVectorLengths )
%     plotframe( __ , Parent=ax )
%     plotframe( __ , Name=Value )
%     hg = plotframe( __ )
% 
%   INPUTS
%   - rotationMatrix        Defines the orientation of the coordinate  
%                           frame, from the origin. 3-by-3 orthogonal 
%                           matrix. Default is zero rotation, i.e., eye(3). 
%   - translationVector     Defines the position of coordinate frame, from 
%                           the origin. 1-by-3 or 3-by-1 numeric vector.
%                           Default is [0 0 0].
%   - basisVectorLengths    Length to plot each arrow (basis) of the
%                           coordinate frame. Scalar, 1-by-3, or 3-by-1 
%                           numeric vector. Default is 1.
%   - Name-Value Arguments
%       + Parent            Axes in which to plot. Scalar axes, group 
%                           (hggroup), or transform (hgtransform) object.
%                           Default is the current axes (gca).
%       + UpdateFrame       With UpdateFrame, the passed plot handles
%                           will be updated with the current parameters, 
%                           rather than creating a new plot. This is more 
%                           efficient and convenient for moving frames.
%                           Handle to an existing frame plot, outputted
%                           from a previous call to plotframe.
%       + MatrixIndexing    Depending on notation, either the columns or 
%                           the rows of rotationMatrix define the 
%                           orientation of the basis vectors. Text scalar, 
%                           either "columnmajor" or "rowmajor". Default is 
%                           row-major.
%       + LabelBasis        Whether the bases should be labelled, e.g., 
%                           "X", "Y", and "Z". Scalar logical. Default is 
%                           false.
%       + Labels            Text with which to label each basis, if 
%                           LabelBasis is enabled. Scalar, 1-by-3, or 
%                           3-by-1 text vector. Default is {'X','Y','Z'}.
%       + BasisColors       Color for each basis vector. Any color format
%                           accepted by MATLAB can be used, e.g., RGB 
%                           triplet [0 0 0], hexadecimal color code 
%                           #000000, or color name 'black' or 'k'. Specify
%                           multiple colors with an M-by-3 matrix where 
%                           each RGB color is a row, or as a 1-by-3 or 
%                           3-by-1 text array. Default is {'r','g','b'}.
%       + TextProperties    Custom properties for the basis labels. 
%                           Name-value arguments stored in a cell array of 
%                           alternating text property names and values, 
%                           e.g., {'FontSize',20,'FontWeight','bold'}.
%       + QuiverProperties  Additional name-value arguments are passed as
%                           properties of the Quiver charts used to plot 
%                           the basis vectors, e.g.,
%                           plotframe( LineStyle="-.", Marker="o" ).
%   OUTPUTS
%   - hg                    Group object (hggroup) containing handles to  
%                           the constituent parts of the coordinate frame 
%                           plot, i.e., the 3 Quiver and optional Text 
%                           objects.
% 
%   EXAMPLES: Please see the file 'examples.mlx' or 'examples.pdf'.
%    
%   Created in 2022b. Compatible with 2020b and later. Compatible with all 
%   platforms. Please cite George Abrahams 
%   https://github.com/WD40andTape/plotframe.
% 
%   See also QUIVER3, QUAT2ROTM, EUL2ROTM, MAKEHGTFORM, PLOTCAMERA
%   Published under MIT License (see LICENSE.txt).
%   Copyright (c) 2023 George Abrahams.
%   - https://github.com/WD40andTape/
%   - https://www.linkedin.com/in/georgeabrahams/
    arguments
        rotationMatrix (3,3) double { mustBeNonNan, mustBeFloat, ...
            mustBeOrthonormal } = eye( 3 )
        translationVector (1,3) double = [ 0, 0, 0 ]
        basisVectorLengths (3,1) double = 1
        Options.MatrixIndexing { mustBeTextScalar, mustBeMember( ...
            Options.MatrixIndexing, [ "columnmajor", "rowmajor"] ) } ...
            = "rowmajor"
        Options.Parent (1,1) { mustBeAxes } = gca
        Options.UpdateFrame matlab.graphics.primitive.Group ...
            { mustBeFrameHG } = matlab.graphics.primitive.Group.empty
        Options.LabelBasis (1,1) logical = false
        Options.Labels (3,1) { mustBeText } = { 'X'; 'Y'; 'Z' }
        Options.BasisColors { mustBeValidColors } = { 'r'; 'g'; 'b' }
        Options.TextProperties cell = {}
        QuiverProperties.?matlab.graphics.chart.primitive.Quiver
        QuiverProperties.AutoScale (1,1) matlab.lang.OnOffSwitchState ...
            { mustBeMember( QuiverProperties.AutoScale, 'off' ) } = 'off'
    end
    isUpdateFrame = ~isempty( Options.UpdateFrame );
    if ~isUpdateFrame
        hg = hggroup( Options.Parent );
    else
        hg = Options.UpdateFrame;
        hg.Parent = Options.Parent;
    end
    [ hQuiver, hText ] = parsehandles( Options.UpdateFrame );
    initgobjects = @( fun ) arrayfun( @(~) fun( 'Parent', hg ), 1:3 );
    
    % If row-major order, the rows of the matrix denote its basis vectors.
    % If column-major order, the columns are its basis vectors.
    if any( isnan( translationVector ) )
        translationVector = [ 0, 0, 0 ];
    end
    if any( isnan( basisVectorLengths ) )
        basisVectorLengths = 1;
    end
    if strncmpi( Options.MatrixIndexing, "columnmajor", 1 )
        rotationMatrix = rotationMatrix';
        Options.MatrixIndexing = "rowmajor";
    end
    basisVectors = rotationMatrix .* basisVectorLengths;
    rgb = validatecolor( Options.BasisColors, 'multiple' );
    quiverProps = namedargs2cell( QuiverProperties );
    if isempty( hQuiver )
        hQuiver = initgobjects( @matlab.graphics.chart.primitive.Quiver );
    end
    set( hQuiver, 'LineWidth', 2, 'MaxHeadSize', 0.4 )
    set( hQuiver, quiverProps{:} )
    set( hQuiver, ...
        { 'XData', 'YData', 'ZData' }, num2cell( translationVector ), ...
        { 'UData', 'VData', 'WData' }, num2cell( basisVectors ), ...
        { 'Color'                   }, num2cell( rgb, 2 ) )
    
    if Options.LabelBasis
        if ~isfield( QuiverProperties, 'Alignment' ) || ...
                strcmpi( QuiverProperties.Alignment, 'tail' )
            textPosition = translationVector + basisVectors;
        elseif strcmpi( QuiverProperties.Alignment, 'center' )
            textPosition = translationVector + basisVectors / 2;
        else % QuiverProperties.Alignment = 'head'
            textPosition = translationVector - basisVectors;
        end
        if isempty( hText )
            hText = initgobjects( @matlab.graphics.primitive.Text );
        end
        set( hText, { 'Position' }, num2cell( textPosition, 2 ), ...
                    { 'String'   }, cellstr( Options.Labels ) )
        if ~isempty( Options.TextProperties )
            try
                set( hText, Options.TextProperties{:} );
            catch ME
                id = "plotframe:InvalidTextProperties";
                msg = "One or more properties or values in the " + ...
                    "TextProperties name-value argument are not valid.";
                ME = addCause( ME, MException( id, msg ) );
                throw( ME )
            end
        end
    elseif ~isempty( hText )
        delete( hText )
    end
    if isUpdateFrame
        drawnow
    end
    if nargout > 0
        varargout = { hg };
    end
end
%% Utility functions
function [ hQuiver, hText ] = parsehandles( hg )
    handlesofclass = @(x) findobj( hg, '-isa', x );
    hQuiver = handlesofclass( 'matlab.graphics.chart.primitive.Quiver' );
    hText = handlesofclass( 'matlab.graphics.primitive.Text' );
end
%% Validation functions
function mustBeOrthonormal( matrix )
    tolerance = 1e-4;
    mxmT = pagemtimes( matrix, "none", matrix, "transpose" );
    eyeDiff = abs( mxmT - eye( 3 ) );
    isOrthonormal = squeeze( all( eyeDiff < tolerance, [1 2] ) );
    if ~isOrthonormal
        id = "plotframe:Validators:MatrixNotOrthonormal";
        msg = sprintf( ...
            "Must be orthonormal (within tolerance %s), i.e., " + ...
            "the basis vectors must be perpendicular and unit length.", ...
            string( tolerance ) );
        throwAsCaller( MException( id, msg ) )
    end
end
function mustBeAxes( x )
    % Validate that x is a valid graphics objects parent, i.e., an axes, 
    % group (hggroup), or transform (hgtransform) object, and has not been 
    % deleted (closed, cleared, etc).
    isAxes = isgraphics( x, "axes" ) || isgraphics( x, "hggroup" ) || ...
         isgraphics( x, "hgtransform" );
    if ~isAxes
        id = "plotframe:Validators:InvalidAxesHandle";
        msg = "Must be handle to graphics objects " + ...
            "parents which have not been deleted.";
        throwAsCaller( MException( id, msg ) )
    end
end
function mustBeFrameHG( hg )
    [ hQuiver, hText ] = parsehandles( hg );
    isFrame = numel( hQuiver ) == 3 && ...
        ( numel( hText ) == 0 || numel( hText ) == 3 );
    if ~isFrame && ~isempty( hg )
        id = "plotframe:Validators:CantUpdateFrame";
        msg = "Must be either an empty Group object or a Group " + ...
            "object returned by a previous call to the function.";
        throwAsCaller( MException( id, msg ) )
    end
end
function mustBeValidColors( colors )
    % validatecolors will either throw an error, which will be handled by
    % the PLOTFRAME arguments block, or return the equivalent RGB colors.
    rgb = validatecolor( colors, 'multiple' );
    nColors = size( rgb, 1 );
    if nColors ~= 1 && nColors ~= 3
        id = "plotframe:Validators:WrongNumberOfColors";
        msg = "Must contain either 1 color (for all 3 axes) or 3 " + ...
            "colors (one for each of the 3 axes).";
        throwAsCaller( MException( id, msg ) )
    end
end