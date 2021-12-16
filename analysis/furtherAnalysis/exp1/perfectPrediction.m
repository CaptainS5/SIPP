% plot the "perfect" prediction of vector average and contrast
clear all; clc; close all
apertureAngles = [-9:3:9]';
apertureSpeed = 10;
dotAngles = [-90, 90];
dotSpeed = 5;

colorPlot = [1, 0, 0; 0, 0, 1];

vecAperture = [cos(apertureAngles/180*pi)*apertureSpeed, sin(apertureAngles/180*pi)*apertureSpeed];
for dotN = 1:length(dotAngles)
    vecDot = [cos(dotAngles(dotN)/180*pi)*dotSpeed, sin(dotAngles(dotN)/180*pi)*dotSpeed];
    vecSum = vecAperture + repmat(vecDot, length(apertureAngles), 1);
    vecContrast = vecAperture - repmat(vecDot, length(apertureAngles), 1);
    dirSum(:, dotN) = atan2(vecSum(:, 2), vecSum(:, 1))/pi*180;
    dirContrast(:, dotN) = atan2(vecContrast(:, 2), vecContrast(:, 1))/pi*180;
end

figure
hold on
plot(apertureAngles, dirSum(:, 1)-apertureAngles, '-', 'color', colorPlot(1, :), 'lineWidth', 2)
plot(apertureAngles, dirSum(:, 2)-apertureAngles, '-', 'color', colorPlot(2, :), 'lineWidth', 2)
xlabel('Aperture angle (deg)')
ylabel('Bias in vector sum')
legend({'dot down', 'dot up'}, 'location', 'best')
saveas(gcf, 'predictSum.pdf')

figure
hold on
plot(apertureAngles, dirContrast(:, 1)-apertureAngles, '-', 'color', colorPlot(1, :), 'lineWidth', 2)
plot(apertureAngles, dirContrast(:, 2)-apertureAngles, '-', 'color', colorPlot(2, :), 'lineWidth', 2)
xlabel('Aperture angle (deg)')
ylabel('Bias in vector difference')
legend({'dot down', 'dot up'}, 'location', 'best')
saveas(gcf, 'predictContrast.pdf')