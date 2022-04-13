# Arduino Magnetometer Calibration Library

Calibrate raw magnetormeter data using matrix inverse or eigenvalue decomposition elliptical fit methods. By the fit algorithm,
hard iron and soft iron offsets are calculated and applied to the
raw magnetometer data. 

<p align="center">
  <img src="https://raw.githubusercontent.com/bunalti/EmbeddedMagnetometerCalibration/main/img/nodistortions.png" width="60%"></a>
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/bunalti/EmbeddedMagnetometerCalibration/main/img/harddistortion.png" width="60%"></a>
</p>

<p align="center">
  <img src="https://raw.githubusercontent.com/bunalti/EmbeddedMagnetometerCalibration/main/img/softdistortions.png" width="60%"></a>
</p>



# Reference

This library is ported from Paul Stoffregen's work. It is simplified to the extent where it only calibrates the magnetometer data.


