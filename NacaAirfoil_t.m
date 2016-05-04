% Andres Vargas
%

classdef NacaAirfoil_t < handle
    %NACAAIRFOIL_T Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        Name;
        Geometry = struct('x',[],'Y',[]);
    end
    
    methods
        function this = NacaAirfoil_t(nacaDigits, chordLength, numberOfPoints)
            
            this.Name = sprintf('NACA%d',nacaDigits);
            
            t = mod(nacaDigits, floor(nacaDigits/100)*100)/100; % the last two digits over 100
            m = floor(nacaDigits/1000)/100; % the first digit over 100
            p =  mod(floor(nacaDigits/100), m*1000)/10;% the second digit over 10
            
            x = linspace(0, chordLength, numberOfPoints/2);
            yt = 5*t*chordLength*(0.2969*sqrt(x/chordLength)+ ...
                (-0.1260)*x/chordLength+(-0.3516)*(x/chordLength).^2+ ...
                0.2843*(x/chordLength).^3+(-0.1015)*(x/chordLength).^4);
                                    
            x1 = x(x<=p*chordLength);
            x2 = x(x>p*chordLength);
            yc = [m*x1/p^2.*(2*p-x1/chordLength) ...
                m*(chordLength-x2)/(1-p)^2.*(1+x2/chordLength-2*p)];
            
            dycdx = [2*m/p^2*(p-x1/chordLength) ...
                2*m/(1-p)^2*(p-x2/chordLength)];
            theta = atan(dycdx);
            
            XU = x-yt.*sin(theta);            
            YU = yc+yt.*cos(theta);            
            XL = x+yt.*sin(theta);
            YL = yc-yt.*cos(theta);
            
            XU(end) = chordLength;
            XL(end) = chordLength;
            
            % coerce XL = XU for vertical alignment
            % note that this changes slightly skin area and area inside 
            % the airfoil
            XL = XU; 
            
            X = [XU(end:-1:1) XL(2:end)];
            Y = [YU(end:-1:1) YL(2:end)];
            
            this.Geometry.X = X';
            this.Geometry.Y = Y';            
        end
    end    
end