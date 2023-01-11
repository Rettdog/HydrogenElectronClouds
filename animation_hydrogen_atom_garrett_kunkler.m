%% Set Up Values

%For best viewing the final scatter plot:
%   1. Click View > Camera Toolbar
%   2. Click Orbit View (leftmost option)
%   3. Click Perspective (option third from the right)

close all;

syms theta phi r

%Constants
hbar = 6.582*10^-16 %ev*s

%Enter the electron's quantum numbers
    n=4;
    l=1;
    ml=linspace(-l, l, 2*l+1);

%Choose whether to display the generated histograms
    displayPlots = false;

%Enter number of random values to use to generate the Monte-Carlo distributions
    numGeneration=1000;

randomDetections = rand(1,numGeneration);

%% ----- Temporal -----

Phi = 1/sqrt(2*pi)*exp(1i*ml*phi);
for i=1:length(ml)
    Theta(i) = generateAssociatedLegendrePolymomial(l,ml(i));
end
Rad = generateRadialSolution(n,l);

RadialPDF = Rad*conj(Rad)*r^2;

GeneralPDFArray =  r.^2.*sin(theta).*conj(sum(Phi.*Theta.*Rad)).*sum(Phi.*Theta.*Rad);

GeneralPDF = vpa(GeneralPDFArray(1))

%Convert General PDF into discrete values
    testAxis = [0.01:0.01:10];
    DiscreteRadialPDF = double(subs(RadialPDF, r, testAxis));
%Numerically solve for most likely radius
    [maxProbValue, maxProbIndex] = max(DiscreteRadialPDF);
    maxRadius = testAxis(maxProbIndex)*3;

%Numerically solve for CDF
    lengthR = 25;
    lengthTheta = 25;
    lengthPhi = 25;
    [rTestAxis, thetaTestAxis, phiTestAxis] = meshgrid(linspace(1, maxRadius, lengthR), linspace(0.01, pi-0.01, lengthTheta), linspace(0.01, 2*pi, lengthPhi));
    discreteGeneralPDF = eval(subs(GeneralPDF, {r, theta, phi}, {rTestAxis, thetaTestAxis, phiTestAxis}));
    totalVolume = sum(discreteGeneralPDF)*(maxRadius/lengthR)*(pi/lengthTheta)*(2*pi/lengthPhi)


%GeneralCDF = @(x,y,z, rad, theta, phi) integral3(GeneralPDF(1), 0, x, 0, y, 0, z)/integral3(GeneralPDF(1), 0, maxRadius, 0, pi, 0, 2*pi)
totalVolume = int(int(int(GeneralPDF(1), rad, 0, maxRadius), theta, 0, pi), phi, 0, 2*pi)
GeneralCDF = @(x,y,z, rad, theta, phi) int(int(int(GeneralPDF(1), r, 0, x), theta, 0, y), phi, 0, z)/totalVolume
%convert to int(int(int())) format

GeneralCDF(1,1,1, r, theta, phi)

%ThetaCDF = @(theta, x) int(ThetaPDF, theta, 0, x)/int(ThetaPDF, theta, 0, pi);

%% ----- Phi -----
Phi = 1/sqrt(2*pi)*exp(1i*ml*phi);

%Azimuthal Probability Distribution Function
    PhiPDF = 1/(2*pi);

%Azimuthal Cumulative Distribution Function
    PhiCDF = @(phi) phi/(2*pi);

%Azimuthal Axis
    phiAxis = linspace(0,2*pi,100);

%Generate histogram of Phi values
    phiHistogram = interp1(PhiCDF(phiAxis), phiAxis, randomDetections, 'pchip');

%Display histogram of Phi values
    if displayPlots
    figure(1);
    title('Phi')
    xlim([0 2*pi])
    histogram(phiHistogram, 100);
    end

%% ----- Theta -----
Theta = generateAssociatedLegendrePolymomial(l,ml);

%Polar Probability Distribution Function
    ThetaPDF = Theta*conj(Theta)*sin(theta)*2*pi;

%Polar Cumulative Distribution Function
    ThetaCDF = @(theta, x) int(ThetaPDF, theta, 0, x)/int(ThetaPDF, theta, 0, pi);

%Polar Axis
    thetaAxis = linspace(0,pi,100);

%Generate histogram of Theta values
    thetaHistogram = interp1(double(arrayfun(@(x) ThetaCDF(theta, x), thetaAxis)), thetaAxis, randomDetections, 'pchip');

%Display histogram of Theta values
    if displayPlots
    figure(2);
    title('Theta')
    xlim([0 pi])
    histogram(thetaHistogram, 100);
    end

%% ----- Radius -----
R = generateRadialSolution(n,l);

%Radial Probability Distribution Function
    RadialPDF = R*conj(R)*r^2;

%Find zeros and domain of the Radial PDF

    %Convert Radial PDF into discrete values
        testAxis = [0.01:0.01:10];
        DiscreteRadialPDF = double(subs(RadialPDF, r, testAxis));

    %Numerically solve for most likely radius
        [maxProbValue, maxProbIndex] = max(DiscreteRadialPDF);
        maxProbRadius = testAxis(maxProbIndex);
    
    %Set the upperbound for radius
        if isempty(maxProbRadius)
            maxRadius=1;
        else
            maxRadius=double(max(maxProbRadius*3, 2));
        end

    %Numerically solve for the zeros of the Radial PDF
        RadialPDFZeros = [];
        
        i=maxProbIndex;
        while i>1
            if DiscreteRadialPDF(i)<=DiscreteRadialPDF(i-1) && DiscreteRadialPDF(i)<=DiscreteRadialPDF(i+1)
                RadialPDFZeros = [RadialPDFZeros, testAxis(i)];
            end
            i=i-1;
        end
        
        RadialPDFZeros = sort(RadialPDFZeros);
    
%Radial Cumulative Distribution Function
    RadialCDF = @(r, x) int(RadialPDF, r, 0, x)/int(RadialPDF, r, 0, maxRadius);

%Radial Axis
    radialAxis = unique(linspace(0,maxRadius,75));

%Generate histogram of Radius values    
    radialHistogram = interp1(double(arrayfun(@(x) RadialCDF(r, x), radialAxis)), radialAxis, randomDetections, 'pchip');

%Display histogram of Radius values
    if displayPlots
    figure(3);
    title('Radius')
    xlim([0 maxRadius])
    histogram(radialHistogram, 100);
    end

%% ----- 3D Scatter -----

%Enter number of simulated electron detections
    numPoints = 20000;

%Simulate points from the Monte-Carlo distributions
    phiLength = length(phiHistogram);
    thetaLength = length(thetaHistogram);
    radialLength = length(radialHistogram);
    
    phiValues = phiHistogram(randi(phiLength, numPoints, 1));
    thetaValues = thetaHistogram(randi(thetaLength, numPoints, 1));
    rValues = radialHistogram(randi(radialLength, numPoints, 1));

%Group the different radii sections by color

%Preallocate the color values to zero
    colorValues = zeros(1, length(rValues));

%Loop through the radius by its zeros
    previousRadius = 0;
    group=1;
    for i = 1:length(RadialPDFZeros)
        group=i;
        currentRadius = RadialPDFZeros(group);
        %Set the color of each radius group based on the loop counter
        currentSection = ((previousRadius <= rValues).*(rValues<currentRadius)).*group;
        colorValues=colorValues+currentSection;
        previousRadius = currentRadius;
    end
%Set the color for the last section
    currentSection = (previousRadius<=rValues).*(group+1);
    colorValues=colorValues+currentSection;

%Convert the Spherical points to Cartesian points
    xValues = rValues.*sin(thetaValues).*cos(phiValues);
    yValues = rValues.*sin(thetaValues).*sin(phiValues);
    zValues = rValues.*cos(thetaValues);

%Plot a scatter graph of the simulated electron detections
    figure(4);
    title('3D Scatter of Hydrogen Atom');
    xlim([-maxRadius maxRadius])
    ylim([-maxRadius maxRadius])
    zlim([-maxRadius maxRadius])
    s = scatter3(xValues, yValues, zValues, 'filled', 'SizeData',5, 'MarkerFaceAlpha',0.1);
    hold on;
    pbaspect([1 1 1]);
    s.CData = colorValues;
    colormap jet
    hold off

%% Functions to Generate Solutions

%Function to generate the Associated Legendre Polynomials for the Polar
%Solution
function P_l_ml = generateAssociatedLegendrePolymomial(l, ml)
    syms theta x

    legendreX = legendreP(l, x);
    diffLegendreX = diff(legendreX, abs(ml));
    diffLegendreTheta = subs(diffLegendreX, cos(theta));

    extraPolynomial = sin(theta)^ml;
    coefficient = ((2*l+1)*factorial(l-ml)/(2*factorial(l+ml)))^(1/2);

    P_l_ml = simplify(diffLegendreTheta*extraPolynomial*coefficient);
    
end

%Function to generate the Radial Solution using Laguerre Polynomials
function R_n_l = generateRadialSolution(n,l)

    syms x r;
    expr = exp(-x)*x^(n+l);
    
    if n-l-1 == 0
        f = expr;
    else
        f = diff(expr,n-l-1);
    end
   
    AssociatedLag = (x^(-2*l-1)*exp(x))/factorial(n-l-1)*f;
    
    %Bohr Radius in nm
        a = 0.0529;

    if n-l-1 == 0
        G = 1;
    else
        G = subs(AssociatedLag, (2*r/(n*a)));
    end
    R_n_l = sqrt((2/(n*a))^3*factorial(n-l-1)/(2*n*factorial(n+l)))*exp(-r/(n*a))*(2*r/(n*a))^l*G;

end