# Unscented Kalman Filter
A more advanced and more accurate method compared to previously done EKF is implemented in this UKF project. UKF does the same task as EKF that is to use lidar and radar data (with noise) to predict the position and velocity of a vehicle.

## Build instructions

Assuming you have 'cmake' and 'make' already:

1. Clone this repo.
2. Make a build directory: `mkdir build && cd build`
3. Compile: `cmake .. && make`
4. Run it: `./UnscentedKF` along with the Term 2 simulator.

## Results

| Input      | Ground truth | RMSE   |
| ---------- | ------------ | ------ |
| position x | 0.09         | 0.0640 |
| position y | 0.10         | 0.0837 |
| velocity x | 0.40         | 0.3335 |
| velocity y | 0.30         | 0.2168 |

The result is also less than the ground truth values.

Here are some visualization graphs of the results.

![](img/visual1.png)

![visual2](img/visual2.png)

![visual3](img/visual3.png)