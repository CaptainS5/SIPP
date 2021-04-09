% create demo for RDK
clc; clear; close all

folders = dir;
for nameN = 4:length(folders)
    demoName = folders(nameN).name;
    
    files = dir([demoName, '\frame*.jpg']);
    for ii = 1:length(files)
        temp = imread([demoName '\frame', num2str(ii), '.jpg']);
        %     temp(temp>20) = temp(temp>20)+0.3*(255-temp(temp>20));
        %     temp(temp<=20) = 0;
        frames{ii} = temp;
    end
    % frame2--fixation, 0.3 + 0.3~0.6 s;
    % frame3, 4--gap, 0.3s
    % frame5--start of RDK
    
    % create the video writer with 1 fps
    writerObj = VideoWriter(['demo_' demoName '.avi']);
    writerObj.FrameRate = 85;
    % open the video writer
    open(writerObj);
    % write the frames to the video
    for u=1:length(frames)
        % convert the image to a frame
        frame = im2frame(frames{u});
        %     if u==2 % fixation
        %         frameN = 40;
        %     elseif u==3
        %         frameN = 17;
        %     else
        %         frameN = 1;
        %     end
        %     for v=1:frameN
        writeVideo(writerObj, frame);
        %     end
    end
    % close the writer object
    close(writerObj);
end

% % for procedure plot
% fixImg = imread('fixation.jpg');
% dispImg = imread('display.jpg');
% flashImg = imread('flash.jpg');
% respImg = imread('response.jpg');
%
% figure
% fixImg = fixImg+0.3*(255-fixImg);
% imshow(fixImg)
% imwrite(fixImg, 'fixation2.jpg')
%
% figure
% dispImg = dispImg+0.3*(255-dispImg);
% imshow(dispImg)
% imwrite(dispImg, 'display2.jpg')
%
% figure
% flashImg = flashImg+0.3*(255-flashImg);
% imshow(flashImg)
% imwrite(flashImg, 'flash2.jpg')
%
% figure
% respImg(respImg>20) = respImg(respImg>20)+0.3*(255-respImg(respImg>20));
% respImg(respImg<=20) = 0;
% imshow(respImg)
% imwrite(respImg, 'response2.jpg')