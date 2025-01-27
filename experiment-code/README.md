# Equipment control/spectrometer operation code
The typical sequence for using the spectrometer chip is
1. Align chip and touch down electrical probes.
2. Hook up santec laser as source and use calibration_acquisition.m to acquire calibration dataset. Do this once with the fiber chuck at 0 (TE) and again at 90 degrees (TM).
3. Use calibration_processing.m to extract calibration results from those two datasets and save in a convenient location.
4. Hook up unknown/test sources.
5. Use spectrometer.m to capture and save individual datasets + OSA references. For some datasets where several measurements are taken in a loop, calibration_acquisition.m can be used instead to automate this.
