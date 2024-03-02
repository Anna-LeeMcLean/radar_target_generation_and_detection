# Radar Target Generation And Detection
The aim of this project is to simulate transmitted and received radar signals from a moving target then detect the target within the simulated signals and determine its displacement from the radar using FFT and CA-CFAR.

## Project Layout

![project layout](https://github.com/Anna-LeeMcLean/radar_target_generation_and_detection/assets/60242931/15b9431f-ba12-4200-956b-4f8a387de3d3)

*Figure 1: Project Workflow*

Two Frequency Modulated Continuous Waves (FMCW) were generated using the system specifications. One wave represented the transmitted signal from a simulated radar and the other the return signal after some time delay.

The System Specs are shown below:

|Specification| Value |
|-|-|
|Frequency| 77 GHz|
| Range Resolution | 1 metre |
| Maximum Detectable Range | 200 metres |
| Maxium Detectable Velocity | 100 m/s |
| Target Initial Distance | 110 meters |
| Target Initial Velocity | -20 m/s |

From these specs, the bandwidth, chirp time and chirp slope for the FMCWs were calculated:

| FMCW Spec | Value |
|-|-|
| Bandwidth | 150e6 GHz |
| Chirp Time | 7.33e-6 secs |
| Chirp Slope | 204.55e11 |

## Range Detection
After generating the transmitted and returning FMCWs, the signals were mixed to determine the beat signal. The FFT over the range values gave a peak +1m from the original initial distance of 110m. The FFT over both the range and doppler values showed a velocity value that was +2 m/s from the original initial velocity of -20 m/s.

![FFT_over_range](https://github.com/Anna-LeeMcLean/radar_target_generation_and_detection/assets/60242931/885d1b68-8753-42d9-b2ef-cdf38ebf00e5)

*Figure 2: 1D FFT showing target distance*

![FFT_over_range_and_velocity](https://github.com/Anna-LeeMcLean/radar_target_generation_and_detection/assets/60242931/3a8d3577-90ea-47ae-84d0-6df01749b3b9)

*Figure 3: 2D FFT showing target distance and velocity*

## Noise Suppression
The 2D FFT plot shows numerous noise spikes. To suppress this noise, CA-CFAR was used to create a dynamic noise threshold that varies with the data.

![RDM_after_CFAR_thresholding](https://github.com/Anna-LeeMcLean/radar_target_generation_and_detection/assets/60242931/e2a6562f-7427-43e6-a3ad-6c629bef0fa9)

*Figure 4: 2D FFT showing target distance and velocity after noise suppression*

The threshold is determined by averaging the training cells around every cells in a Range-Doppler Map while avoiding guard cells so that the average doesn't consider the cell under test.
![image](https://github.com/Anna-LeeMcLean/radar_target_generation_and_detection/assets/60242931/a49343fa-dc2f-41fa-bcaa-c1d24cee2d6e)

*Figure 5: CA-CFAR grid used to determine noise threshold for a Cell Under Test (CUT)*

Figure 5 shows the training and gurad cells around the cell under test (CUT). After averaging the training cells, an offset is added to the threshold which allows for added padding and easy tuning of the threshold value.  

| CA-CFAR Specs | Value |
|-|-|
| Number of training cells for range bin of RDM | 10 |
| Number of training cells for doppler bin of RDM | 8 |
| Number of guard cells for range bin of RDM | 4 |
| Number of guard cells for doppler bin of RDM | 4 |
| Offset | 4.5 |

The offset was increased from an initial value of 1.5 to 4.5 after observing that most of the noise spikes were still present. Increasing the threshold by this offset allowed for larger noise spikes to be suppressed without cutting off the real signal values. All signal values at the CUT were suppressed to 0 if they were below the threshold while values which were above the threshold were set to 1. The new values were stored at their same locations in a new matrix which was the same size as the original RDM matrix. While the CFAR grid is at the edges of the Range Doppler map, the cells within the training and guard cell areas cannot be tested as CUTs because they would not have sufficient training/guard cells around them. Initializing the new RDM matrix as zeros allows these cells that are not processed to remain as zeros.
