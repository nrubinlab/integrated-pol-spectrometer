% Prompts UI to select processed calibration data and performs basic checks
function [acq_struct, TE_calib_file, TM_calib_file, amp_calib_file] = prompt_load_calibration(options)
    arguments
        options.calib_dirname = []
    end
    % For simplicity calibration files MUST be named "TE.mat" and "TM.mat" in a
    % single folder.
    if(isempty(options.calib_dirname))
        calib_dirname = uigetdir();
    else
        calib_dirname = options.calib_dirname;
    end
    
    TE_calib_file = fullfile(calib_dirname, 'TE.mat');
    TM_calib_file = fullfile(calib_dirname, 'TM.mat');
    % deprecated amplitude calibration functionality not used in this work.
    amp_calib_file = fullfile(calib_dirname, 'amplitude.mat');
    if(isfile(TE_calib_file))
        load(TE_calib_file, "acq_struct");
        TE_acq_struct = acq_struct;
    else
        error("No 'TE.mat' calib file in selected folder!");
    end
    
    if(isfile(TM_calib_file))
        load(TM_calib_file, "acq_struct");
        if(~isequal(TE_acq_struct, acq_struct))
            error("TE and TM acquisition structures are not identical!");
        end
    else
        error("No 'TM.mat' calib file in selected folder!");
    end

    if(~isfile(amp_calib_file))
        disp("No 'amplitude.mat' calib file in selected folder, returning NaN.");
        amp_calib_file = NaN;
    end
    disp("Calibration import successful.");
end

