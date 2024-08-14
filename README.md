# HL-Detection-using-TSA

## Authors
[Ahmet Agaoglu](https://github.com/Ahmet-Agaoglu), [Nezih Topaloglu](https://github.com/nezihtopaloglu)

## Description
This repository contains MATLAB code for Horizon Line (HL) Detection using dynamic ROI generation based on Time Series Models (ARIMA and GARCH). The approach focuses on fast and accurate HL detection in maritime videos, particularly in dynamic environments.

## Files Included

The repository includes a MATLAB script and two folders:

- **HL_Detect_TSA.m**: This script processes an input video and outputs the estimated Horizon Line state matrix. Each video frame is annotated with a parallelogram ROI (in yellow) and the estimated HL (in red). Run this script to execute the main functionality.

- **VIS_Onboard**: an empty folder for Onboard Singapore Maritime Dataset videos.

- **Buoy**: an empty folder for Buoy Dataset videos.

## Getting Started

### Installation:
1. Download the repository to your system.
2. Download the Singapore Maritime Dataset On-Board videos into `~/VIS_Onboard` from [this link](https://drive.google.com/file/d/0B43_rYxEgelVb2VFaXB4cE56RW8/view?resourcekey=0-67PrivAOYTGyWxAO_-2n1A).
3. Download the Buoy Dataset videos into `~/Buoy` from [this link](https://drive.google.com/file/d/0B43_rYxEgelVVngtMVBpWGFqckE/view?resourcekey=0-zBgpYCkkblxPZocaf8NU5w).
4. Set the downloaded repository folder as the MATLAB working directory. Add the video folders to the MATLAB path.

### Usage

To run the code, call the function `HL_Detect_TSA` with the input video name.

## Citation
If you use this code in your research, please cite the following paper:
```
@article{Agaoglu2024,
  title={Enhancing Maritime Horizon Line Detection Using Time Series Analysis},
  author={Ahmet Agaoglu and Nezih Topaloglu},
  journal={Under Review},
  year={2024}
}
```
## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

