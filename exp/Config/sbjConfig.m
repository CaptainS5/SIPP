function [sbj] = sbjConfig(const)
% Type in the subject information and store in sbj:
%   PK, 21 / 03 / 2019
% modified by XiuyunWu 24/02/2020

if const.startExp ~=0
   sbj.name   = input(sprintf('\n\tID: '),'s');
   sbj.block  = input(sprintf('\n\tFrom Block: '));
   sbj.trial  = input(sprintf('\n\tFrom Trial (within the current block): '));
   if sbj.block == 1 && sbj.trial==1
      sbj.age    = input(sprintf('\n\tAge: '));
      sbj.sex    = input(sprintf('\n\tSex: '),'s');
      sbj.hand   = input(sprintf('\n\tHandedness: '),'s');
   end
else
   sbj.name   = 'Anon';
   sbj.age    = 0;
   sbj.sex    = 'Anon';
   sbj.hand   = 'Anon';
   sbj.block  = 1;
   sbj.trial  = 1;
end
currentTime = clock;
sbj.date = sprintf('%d-%d-%d_%d%d', currentTime(1:5));

% create a folder in data: e.g. data/PK
currentFolder = pwd;
if const.startExp==-1 % practice blocks
    sbj.sbjFolder = ['../data/',sbj.name, '/practice'];
else
    sbj.sbjFolder = ['../data/',sbj.name];
end
if ~exist(sbj.sbjFolder)
    mkdir(sbj.sbjFolder);
end
cd([sbj.sbjFolder, '/'])
sbj.sbjFolder = pwd;

% create subfolders 
% for saving rdk seeds
sbj.rdkFolder = [sbj.sbjFolder '/rdkSeeds'];
if ~exist(sbj.rdkFolder)
    mkdir(sbj.rdkFolder);
end
% % for each block and task:
% sbj.currentFolder = ['data/',sbj.name, '/', sbj.task, '_', sbj.block];
% if ~exist(sbj.currentFolder)
%     mkdir(sbj.currentFolder)
% else
%     overwrite = input(sprintf('\n\tDo you want to delete and overwrite old fodler? (1 - NO, 2 - YES): '));
%     if overwrite == 1
%         sbj.currentFolder = ['data/',sbj.name, '/', sbj.task, '_', sbj.block, '_new'];
%         if exist(sbj.currentFolder)
%             disp('EXP: ERROR - CHECK YOUR FOLDERS!')
%             error;
%         else
%             mkdir(sbj.currentFolder)
%         end
%     elseif overwrite == 2
%         mkdir(sbj.currentFolder)
%     end
% end
cd(currentFolder)
end