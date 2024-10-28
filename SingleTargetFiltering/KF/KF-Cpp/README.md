### Kalman Filter (KF)
 Kalman filter implementation using C++ with Boost based on the following
 introductory paper:
     http://www.cs.unc.edu/~welch/media/pdf/kalman_intro.pdf

It predicts an object's location in 2D space i.e. it takes positions (xy coordinates) 
as input and returns estimated positions and speed (xy coordinates and v_x v_y speed). 
It uses a constant velocity (CV) motion model. 

### Dependencies
- C++11
- boost ublas

### Running Example Code
```bash
mkdir build
cd build
cmake ..
make
./KalmanFilter_demo
```

### Disply KF results
It uses Python to display the ground truth and estimated positions after saving them in a file.
After running the ./KF_demo above, use the following to display its result.
```bash
python3 KalmanFilter_display.py
```



