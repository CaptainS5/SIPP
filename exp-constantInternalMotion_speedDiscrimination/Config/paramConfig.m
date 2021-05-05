function [trialCondition, trialData] = paramConfig(const,sbj)
% =========================================================================
% paramConfig(const,sbj)
% =========================================================================
% Function created to set up all Experimental Parameters (constants and
% randomized variables). Performs permutation for all randomized variables.
%
% -------------------------------------------------------------------------
% Input:
% const:     structure containing different constant settings
% sbj:       structure containing subject information
% -------------------------------------------------------------------------
% Output:
% trialCondition: structure containing all experimental parameters for all trials.
%   should be used as a look-up table during the experiment
% -------------------------------------------------------------------------
%%
if sbj.block == 1 && sbj.trial==1 % ~exist([sbj.sbjFolder ,'/trialData.mat'],'file')
    rng('shuffle') % just do it once per experiment
    
    trialCondition = table;
    % (1) SET UP MEASURE VARIABLES                                              % These will be filled in during the experiment
    % now just write into trialData directly during the experiment
    trialCondition.tMainSync          = zeros(const.numTrials, 1);                   % Time (GetSecs) at trial start, also fixation on.
    trialCondition.tRDKon        = NaN(const.numTrials, 1);
    trialCondition.tRDKoff  = NaN(const.numTrials, 1);                     % actual measured time that the target appeared/disappeared
    %     trialCondition.tResponse          = NaN(const.numTrials, 1);               % also the end of the trial
    %     %     trialData.t_start_VBL        = NaN(const.numTrials, 5);
    %     %     trialData.t_rdkOn_VBL         = NaN(const.numTrials, 5);
    %     %     trialData.t_rdkOff_VBL        = NaN(const.numTrials, 5);
    %     %     trialData.t_response_VBL        = NaN(const.numTrials, 5);
    %     trialCondition.choice          = NaN(const.numTrials, 1); % -1=up, 1=down
    % Sync Times Setup:
    
    
    % (2) RANDOMIZED VARIABLES
    % randVar.continuous: pick from the given range uniformly for each trial;
    % such as fixation duration
    % randVar.trial: condition in each trial, discrete numbers
    % randVar.block: blocked condition, randomize the order of all blocks, discrete
    % numbers
    % randVar.blockPartial: blocked condition, not completely randomized,
    % but present blocks with the same conditions together; only randomize
    % the order of all combination of conditions
    % if both blockPartial and block exist, "block" should be randomized
    % within "blockPartial"; all blockPartial will have the same "block"
    % conditions in total
    randVar.continuous.fixationDuration = [1 1.5]; % Intial fixation time (random between 300 and 700 ms)
    randVar.continuous.rdkDuration = const.rdk.duration;
    %     randVar.trial.dotDirSD = const.rdk.dotDirSD;
    %     randVar.trial.rdkCoh = const.rdk.coh;
    randVar.trial.rdkApertureDir = const.rdk.apertureDir;
    randVar.trial.rdkApertureSpeed = const.rdk.apertureSpeed;
    randVar.trial.rdkInternalCons = const.rdk.internalCons;
    %     randVar.trial.rdkInternalDir = const.rdk.internalDir; % upwards is minus, and downwards is plus
    
    %     if const.startExp~=-1 % make accurate blocks first, and then fast, just so that it's easier to learn
    %         randVar.blockPartial.instruction = [0 1]; % 0-fast, 1-accurate;
    %         randVar.blockPartial.eyeCondition = [0 1]; % 0-fixation, 1-pursuit
    %         % if simply want to randomly interleave all blocks,
    %         % comment the above lines under "if" and uncomment below
    %         % % ============================================================
    %         % randVar.block.instruction = [0 1]; % 0-fast, 1-accurate; if no blocked conditions, set to []
    %         % randVar.block.eyeCondition = [0 1]; % 0-fixation, 1-pursuit
    %         % % ============================================================
    %     end
    
    if const.startExp>=0
        % perform randomization:
        trialCondition = [trialCondition pick_paramsAll(randVar, const.numTrialsPerBlock)];
        
        % (3) CONSTANT VARIABLES                                                    % Set some constants that are used similarly in all blocks
        % insert standard trials every ten trials (trial 1, 11, ...)
        standardN = size(trialCondition, 1)/const.standardInterval;
        for insertN = 1:standardN
            % set the parameters for the standard trial to be inserted
            tempCon = table;
            tempCon.tMainSync = 0;                   % Time (GetSecs) at trial start, also fixation on.
            tempCon.tRDKon = NaN;
            tempCon.tRDKoff = NaN;
            tempCon.trialConditionIdx = const.standardInterval*(insertN-1)+insertN; % need to also add the interted trials before...
            % need to add the later trial idx in trialCondition by one
            idxTemp = find(trialCondition.trialConditionIdx==tempCon.trialConditionIdx);
            trialCondition.trialConditionIdx(idxTemp:end, 1) = trialCondition.trialConditionIdx(idxTemp:end, 1)+1;
            
            tempCon.blockN = trialCondition.blockN(idxTemp, 1); % still within the same block
            tempCon.rdkApertureDir = 0;
            tempCon.rdkApertureSpeed = 10; % the standard speed
            tempCon.rdkInternalCons = -1; % static pattern
            tempCon.fixationDuration = max(randVar.continuous.fixationDuration); %min(randVar.continuous.fixationDuration)+rand*(max(randVar.continuous.fixationDuration)-min(randVar.continuous.fixationDuration));
            tempCon.rdkDuration = min(randVar.continuous.rdkDuration)+rand*(max(randVar.continuous.rdkDuration)-min(randVar.continuous.rdkDuration));
            % now inset into trialCondition
            if idxTemp==1
                trialCondition = [tempCon; trialCondition(idxTemp:end, :)];
            else
                trialCondition = [trialCondition(1:idxTemp-1, :); tempCon; trialCondition(idxTemp:end, :)];
            end
        end
        
    else % practice block, basically only standard trials
        trialCondition.trialConditionIdx(:, 1) = 1:size(trialCondition, 1)';
        trialCondition.blockN(:, 1) = 1; % still within the same block
        trialCondition.rdkApertureDir(:, 1) = 0;
        trialCondition.rdkApertureSpeed(1:10, 1) = 10; % the standard speed
        trialCondition.rdkInternalCons(1:10, 1) = -1; % static pattern
        trialCondition.fixationDuration(:, 1) = min(randVar.continuous.fixationDuration) + ...
            rand(size(trialCondition, 1), 1)*(max(randVar.continuous.fixationDuration) - min(randVar.continuous.fixationDuration));
        trialCondition.rdkDuration(:, 1) = min(randVar.continuous.rdkDuration) + ...
            rand*(max(randVar.continuous.rdkDuration) - min(randVar.continuous.rdkDuration));
        if size(trialCondition, 1)>10 % random test conditions
            for ii = 1:size(trialCondition, 1)-10
                trialCondition.rdkApertureSpeed(10+ii, 1) = const.rdk.apertureSpeed(randi(length(const.rdk.apertureSpeed))); % the standard speed
                trialCondition.rdkInternalCons(10+ii, 1) = const.rdk.internalCons(randi(length(const.rdk.internalCons)));
            end
        end
    end
    trialCondition.rdkDistance = trialCondition.rdkApertureSpeed.*trialCondition.rdkDuration;
    
    % (4) SAVE TRIALDATA
    save([sbj.sbjFolder ,'/trialCondition.mat'], 'trialCondition');
    trialData = table; % simply initialize trialData for later use
else
    load([sbj.sbjFolder ,'/trialCondition.mat']);
    load([sbj.sbjFolder ,'/trialData.mat']) % trialData will be generated once the experiment started; load this as well
end

