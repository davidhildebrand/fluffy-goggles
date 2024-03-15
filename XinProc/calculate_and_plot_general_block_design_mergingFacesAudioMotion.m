% You first need to run the general calculate_and_plot_face_contrast with
% some face-trial-data!!! 

fps = 5;
saving = false;
time_delay_hemo = 4;

n_trials = 22;
time_stim_start = 4;
time_stim_end = 14;

test_additional_hemo_delays = false;

[filenames, folder, ~] = uigetfile( 'D:\XINTRINSIC\*_P1.mat', 'Select MOTION/DOTS/MT preprocessed files', 'MultiSelect', 'On');

for i_extra_delay = 0 %1:15 % Only loops through the 15 if test_additional_hemo_delays = true

    result_map_files = cell(length(filenames),1);
    for i_file = 1:length(filenames)
        disp(filenames{i_file})
        load(strcat(folder,filenames{i_file}))

        P.ProcDataMat = squeeze(P.ProcDataMat);

        dims = size(P.ProcDataMat);
        total_frames = dims(1) * dims(4);
        frames_per_trial = total_frames / n_trials;

        concatenated_trials = zeros(dims(2),dims(3),total_frames);
        for i = 1:dims(1)
            this_trial_frame_start = (i-1)*dims(4)+1;
            this_trial_frame_end = i*dims(4);
            concatenated_trials(:,:,this_trial_frame_start:this_trial_frame_end) = P.ProcDataMat(i,:,:,:);
        end

        resorted_trials = zeros(dims(2), dims(3), n_trials, frames_per_trial);
        for i = 1:n_trials
            this_trial_frame_start = (i-1)*frames_per_trial+1;
            this_trial_frame_end = i*frames_per_trial;
            resorted_trials(:,:,i,:) = concatenated_trials(:,:,this_trial_frame_start:this_trial_frame_end);
        end

        resorted_trials = 100*resorted_trials./max(max(max(max(resorted_trials))));

        average_trial = squeeze(mean(resorted_trials, 3));

        frame_response_start = fps * (time_stim_start  + time_delay_hemo);
        frame_response_end   = fps * (time_stim_end    + time_delay_hemo);
        response_bool = false(size(average_trial,3),1);
        response_bool(frame_response_start:frame_response_end) = true;
        response = mean(average_trial(:,:, response_bool),3);
        baseline = mean(average_trial(:,:,~response_bool),3);

        baseline_minus_response = baseline - response;

        figure;
        imagesc(baseline_minus_response);
        %imagesc(baseline_minus_response,[0.1,0.75]); % through-window
        %imagesc(baseline_minus_response,[0.0,0.1]); % through-skull
        title(strcat(filenames{i_file}, ' - ', num2str(n_trials), ' cycles'))
        axis equal; colorbar;
        axis off;

        result_map_files(i_file) = {baseline_minus_response};
    end

    concatenated_maps = cat(3, result_map_files{:});
    average_map = mean(concatenated_maps, 3);

    figure;
    fm = imagesc(average_map,[0.0,0.6]); % through-window
    %imagesc(average_map,[0.0,0.4]); % through-skull
    axis equal; colorbar; axis off;
    %title(strcat('Merged map:', num2str(n_trials * length(filenames)), ' cycles'))
    title('resp amp diff, dots moving vs static averaged over all reps');
    if saving
        monkey = 'Scrooge';
        datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        saveas(fm,[save_path filesep datestrf '_' monkey '_AmpMeanDiffDots_' num2str(cycle_num) 'reps_FigSaveAs.png'])
        saveas(fm,[save_path filesep datestrf '_' monkey '_AmpMeanDiffDots_' num2str(cycle_num) 'reps_FigSaveAs.svg'])
    end
    
    if ~test_additional_hemo_delays
        break
    else
        average_for_this_hemo_delay(i_extra_delay) = mean(mean(average_map))
        time_delay_hemo = time_delay_hemo + 0.2;
    end
    
end

average_map_thresholded_m = average_map;

n_trials = 20;
time_stim_start = 5;
time_stim_end = 14.8;

test_additional_hemo_delays = false;

[filenames, folder, ~] = uigetfile( 'D:\XINTRINSIC\*_P1.mat', 'Select AUDIO/SOUND preprocessed files', 'MultiSelect', 'On');

for i_extra_delay = 0 %1:15 % Only loops through the 15 if test_additional_hemo_delays = true

    result_map_files = cell(length(filenames),1);
    for i_file = 1:length(filenames)
        disp(filenames{i_file})
        load(strcat(folder,filenames{i_file}))

        P.ProcDataMat = squeeze(P.ProcDataMat);

        dims = size(P.ProcDataMat);
        total_frames = dims(1) * dims(4);
        frames_per_trial = total_frames / n_trials;

        concatenated_trials = zeros(dims(2),dims(3),total_frames);
        for i = 1:dims(1)
            this_trial_frame_start = (i-1)*dims(4)+1;
            this_trial_frame_end = i*dims(4);
            concatenated_trials(:,:,this_trial_frame_start:this_trial_frame_end) = P.ProcDataMat(i,:,:,:);
        end

        resorted_trials = zeros(dims(2), dims(3), n_trials, frames_per_trial);
        for i = 1:n_trials
            this_trial_frame_start = (i-1)*frames_per_trial+1;
            this_trial_frame_end = i*frames_per_trial;
            resorted_trials(:,:,i,:) = concatenated_trials(:,:,this_trial_frame_start:this_trial_frame_end);
        end

        resorted_trials = 100*resorted_trials./max(max(max(max(resorted_trials))));

        average_trial = squeeze(mean(resorted_trials, 3));

        frame_response_start = fps * (time_stim_start  + time_delay_hemo);
        frame_response_end   = fps * (time_stim_end    + time_delay_hemo);
        response_bool = false(size(average_trial,3),1);
        response_bool(frame_response_start:frame_response_end) = true;
        response = mean(average_trial(:,:, response_bool),3);
        baseline = mean(average_trial(:,:,~response_bool),3);

        baseline_minus_response = baseline - response;

        figure;
        imagesc(baseline_minus_response);
        %imagesc(baseline_minus_response,[0.1,0.75]); % through-window
        %imagesc(baseline_minus_response,[0.0,0.1]); % through-skull
        title(strcat(filenames{i_file}, ' - ', num2str(n_trials), ' cycles'))
        axis equal; colorbar;
        axis off;

        result_map_files(i_file) = {baseline_minus_response};
    end

    concatenated_maps = cat(3, result_map_files{:});
    average_map = mean(concatenated_maps, 3);

    figure;
    fm = imagesc(average_map,[0.0,0.6]); % through-window
    %imagesc(average_map,[0.0,0.4]); % through-skull
    axis equal; colorbar; axis off;
    %title(strcat('Merged map:', num2str(n_trials * length(filenames)), ' cycles'))
    title('resp amp diff, dots moving vs static averaged over all reps');
    if saving
        monkey = 'Scrooge';
        datestrf = [datestr(now,'yyyymmdd') 'd' datestr(now,'HHMMSS') 't'];
        saveas(fm,[save_path filesep datestrf '_' monkey '_AmpMeanDiffDots_' num2str(cycle_num) 'reps_FigSaveAs.png'])
        saveas(fm,[save_path filesep datestrf '_' monkey '_AmpMeanDiffDots_' num2str(cycle_num) 'reps_FigSaveAs.svg'])
    end
    
    if ~test_additional_hemo_delays
        break
    else
        average_for_this_hemo_delay(i_extra_delay) = mean(mean(average_map))
        time_delay_hemo = time_delay_hemo + 0.2;
    end
    
end


average_map_thresholded_a = average_map;




average_map_thresholded_f = D.R.AmpMeanDiff_tmp;

blurred = true
if blurred
    average_map_thresholded_m2 = imgaussfilt(average_map_thresholded_m,2);
    average_map_thresholded_a2 = imgaussfilt(average_map_thresholded_a,1.5);
    average_map_thresholded_f2 = imgaussfilt(average_map_thresholded_f,2);
else
    average_map_thresholded_m2 = average_map_thresholded_m;
    average_map_thresholded_a2 = average_map_thresholded_a;
    average_map_thresholded_f2 = average_map_thresholded_f;
end


%%% Window
bgim = gray2rgb(R.windowavg*0.000004);
fgim = zeros(size(bgim));

% % Crumpet
% average_map_thresholded_m2(average_map_thresholded_m2 < 0.40) = 0;
% average_map_thresholded_a2(average_map_thresholded_a2 < 0.50) = 0;
% average_map_thresholded_f2(average_map_thresholded_f2 < 0.0065) = 0;

% Curly
average_map_thresholded_m2(average_map_thresholded_m2 < 0.65) = 0;
average_map_thresholded_a2(average_map_thresholded_a2 < 0.50) = 0;
average_map_thresholded_f2(average_map_thresholded_f2 < 0.001) = 0;

%fgim(:,:,1) = squeeze(average_map_thresholded_f2)*100000;
%fgim(:,:,2) = squeeze(average_map_thresholded_m2)*100000;
fgim(:,:,3) = squeeze(average_map_thresholded_a2)*100000;


fgim(bgim < 0.6) = 0; % This is for removing noise from the edges of the chamber/very dark areas

% %% Skull
% bgim = gray2rgb(R.windowavg*0.000008);
% fgim = zeros(size(bgim));
% 
% % Curly
% average_map_thresholded_m2(average_map_thresholded_m2 < 0.27) = 0;
% average_map_thresholded_a2(average_map_thresholded_a2 < 0.037) = 0;
% average_map_thresholded_f2(average_map_thresholded_f2 < 0.00025) = 0;
% 
% fgim(:,:,1) = squeeze(average_map_thresholded_f2)*100000000;
% fgim(:,:,2) = squeeze(average_map_thresholded_m2)*500000000;
% fgim(:,:,3) = squeeze(average_map_thresholded_a2)*100000000;
% 
% fgim(bgim < 1.5) = 0; % This is for removing noise from the edges of the chamber/very dark areas
% 



blended = imblend(fgim, bgim, 0.8, 'normal');
figure
f1h1 = imshow(blended(:,80:400,:));
axis equal; axis off;





