% FUNCTION to update trial information in a text box when clicking through
% trials using viewEyeData.m
% history
% 07-2012       JE created updatePlots.m
% 2012-2018     JF added stuff to and edited updatePlots.m
% 16-07-2018    JF commented to make the script more accecable for future
%               VPOM students
% for questions email jolande.fooken@rwth-aachen.de
%
% input: trial --> structure containing relevant current trial information
%                  (remember by now we have stored saccade information in
%                   trial as well)
% output: this will plot in the open figure

function [] = updateText(trial, fig)

screenSize = get(0,'ScreenSize');
name = evalin('base', 'name');
% chose position of the text box
xPosition = 10; %screenSize(3)*2/3-50;
yPosition = 270; %screenSize(4)-screenSize(4)/2; %screenSize(4)*2/5-50;
% how large should it be?
verticalDistance = 20;
width = 100;
height = 20;
textblock = 0;

% now we can add all infos we want in a text box each
subjectIdText = uicontrol(fig,'Style','text',...
    'String', ['subject ID: ' trial.log.subject],...
    'Position',[xPosition yPosition width height],...
    'HorizontalAlignment','left'); %#ok<*NASGU>

textblock = textblock+1;
trialNoText = uicontrol(fig,'Style','text',...
    'String', ['Trial ' num2str(trial.log.trialNumber), ', block ', num2str(trial.log.blockN)],...
    'Position',[xPosition yPosition-textblock*verticalDistance width height],...
    'HorizontalAlignment','left');

% if trial.log.eyeCondition==0
%     eyeConditionName = 'Fixation';
% else
%     eyeConditionName = 'Pursuit';
% end
% if trial.log.instruction==0
%     instructionName = 'Fast';
% else
%     instructionName = 'Accurate';
% end
% textblock = textblock+1;
% trialNoText = uicontrol(fig,'Style','text',...
%     'String', [eyeConditionName, ', ', instructionName],...
%     'Position',[xPosition yPosition-textblock*verticalDistance width height],...
%     'HorizontalAlignment','left');

if trial.log.rdkApertureDir==0
    rdkDir = 'right';
%     angle = trial.log.rdkApertureAngle;
else
    rdkDir = 'left';
%     angle = 180-trial.log.rdkApertureAngle;
end
% if trial.log.rdkInternalDir>0
%     sdText = 'up';
% else
%     sdText = 'down';
% end
textblock = textblock+1;
trialNoText = uicontrol(fig,'Style','text',...
    'String', ['Aperture: ', rdkDir, ', ', num2str(trial.log.rdkApertureSpeed)],...
    'Position',[xPosition yPosition-textblock*verticalDistance width height],...
    'HorizontalAlignment','left');

if trial.log.rdkInternalDir==0
    sdText = 'right';
else
    sdText = 'left';
end
textblock = textblock+1;
trialNoText = uicontrol(fig,'Style','text',...
    'String', ['Internal: coh ', num2str(trial.log.rdkCoh), ',', sdText],...
    'Position',[xPosition yPosition-textblock*verticalDistance width height],...
    'HorizontalAlignment','left');

if trial.log.response>0
    choice = 'faster';
else
    choice = 'slower';
end
textblock = textblock+1;
trialNoText = uicontrol(fig,'Style','text',...
    'String', ['Response: ', choice],...
    'Position',[xPosition yPosition-textblock*verticalDistance width height],...
    'HorizontalAlignment','left');

end


