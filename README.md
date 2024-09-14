# Hybrid-TDOA-Multi-Calib
In this paper, to enhance calibration accuracy, we propose to incorporate the time difference of arrival measurements between adjacent sound events (TDOA-S) with respect to the microphone arrays. More specifically, we propose a two-stage approach including an initial value estimation (IVE) procedure and the final joint optimization step. The VLE stage initializes all parameters except for microphone orientations using hybrid TDOA (i.e., TDOA-M and TDOA-S), odometer data from a moving robot carrying a speaker, and DOA and refines the microphone orientations. The final joint optimization step estimates multiple microphone array locations, orientations, time offsets, clock drift rates, and sound source locations simultaneously. Both simulation and experiment results show that for scenarios with moderate TDOA noise levels, our approach outperforms existing methods in terms of accuracy.

## Calibration Scenario
<img src="https://github.com/Chen-Jacker/Hybrid-TDOA-Calib/blob/main/calibration_scenario.gif" width="500px">
The full video with high quality is avaliable at https://youtu.be/pmVgjJdjW2E.

## Audio Accessment
Link in OneDrive: https://1drv.ms/f/s!AilTdY3K-LzJgbdgoZZSl9-d8883ow?e=gFOias. After download it, place folder named "Audio" at the same level as folder "displacement". If the link has expired, send me an email to remind me.

## How to implement simulations
- `sim_main.m` achieves simulations in paper.
  1. **Simulations Set Up** generate data of simulation scenario under three levels of TDOA noises, DOA noises, three sound source trajectories. In each TDOA noise and DOA noise and specific trajectory, Monte Carlo performs `nums` times. If you wanna **Quick Simulations**, set `nums` to smaller, such as five.
  2. After running first step, you can run **Zhang's Simulation** in `sim_main.m` and Wang's approach in file named `Multi_Mic_Arrays_Calibration.py`. Note that declare `category='Simulation'`.
  3. Finally, just run the bottom part in `sim_main.m`. Then, IVE and final results of our method and wang's method will be shown.
- `real_main.m` achieves real-world experiments in paper.
  1. Download audio shared in this link and put audio into folder named `audio`.
  2. The first part and second part in `real_main.m` generate TDOA-S and TDOA-M measurements, and real-world experiment data respectively. run `SRP-PHAT-DOA.py` to compute DOA.
  3. Subsquencely, run the third part in `real_main.m` that is my approach and Wang's approach in file named `Multi_Mic_Arrays_Calibration.py`. Note that declare `category='Real_world'`.
  4. Finally, run the final part of `real_main.m` to output the IVE and final estimation results.
## Details
- My email: 12332644@mail.sustech.edu.cn.
