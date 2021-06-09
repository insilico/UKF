# Unscented Kalman Filter (UKF)

### Examples
`inst/model_voxels_2param.R` example script for applying UKF to a pair of voxel time series in `data/voxel_A&B_data.txt` and optimizing coupling parameters for an ode model of synchronization. todo: bind data to the library as rda?

`inst/test_simulate.R` example script to simulate pair of voxel data. todo: fix runge-kutta.
### Abstract

#### Websites

#### Related References

### To install:

    >library(devtools)
    >install_github("insilico/UKF")  
    >library(UFK)

### Dependencies
Used for some pre-process smoothing of time series, but not required for UKF. 
```
install.packages(c('KernSmooth'))
```
#### Contact
[brett.mckinney@gmail.com](brett.mckinney@gmail.com)
