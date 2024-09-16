# HL-Detection-using-TSA [![DOI](https://zenodo.org/badge/842448544.svg)](https://zenodo.org/doi/10.5281/zenodo.13768683)


## Authors
[Ahmet Agaoglu](https://github.com/Ahmet-Agaoglu), [Nezih Topaloglu](https://github.com/nezihtopaloglu)

## Description
This repository contains MATLAB code for Horizon Line (HL) Detection using dynamic ROI generation based on Time Series Models (ARIMA and GARCH). The approach focuses on fast and accurate HL detection in maritime videos, particularly in dynamic environments.

## Files Included
The repository includes a MATLAB script and two folders:

- **HL_Detect_TSA.m**: This script processes an input video and outputs the estimated Horizon Line state matrix. Each video frame is annotated with a parallelogram ROI (in yellow) and the estimated HL (in red). Run this script to execute the main functionality.
- **VIS_Onboard**: An empty folder for Onboard Singapore Maritime Dataset videos.
- **Buoy**: An empty folder for Buoy Dataset videos.

## Datasets
The code requires the following datasets:

- **Singapore Maritime Dataset (Onboard Segment)**: Download the Onboard videos into the `~/VIS_Onboard` folder from [this link](https://drive.google.com/file/d/0B43_rYxEgelVb2VFaXB4cE56RW8/view?resourcekey=0-67PrivAOYTGyWxAO_-2n1A). If using this dataset, please cite: D.K. Prasad, D. Rajan, L. Rachmawati, E. Rajabaly, and C. Quek, "Video Processing from Electro-optical Sensors for Object Detection and Tracking in Maritime Environment: A Survey," IEEE Transactions on Intelligent Transportation Systems (IEEE), 18(8), 1993 - 2016, 2017.
- **Buoy Dataset**: Download the videos into the `~/Buoy` folder from [this link](https://drive.google.com/file/d/0B43_rYxEgelVVngtMVBpWGFqckE/view?resourcekey=0-zBgpYCkkblxPZocaf8NU5w). If using this dataset, please cite: S. Fefilatyev, D. Goldgof, M. Shreve, and C. Lembke, “Detection and tracking of ships in open sea with rapidly moving buoy-mounted camera system,” Ocean Eng. 54, 1–12 (2012).

You can add additional datasets to the respective folders as needed.

## Installation
1. Download the repository to your system.
2. Download the Singapore Maritime Dataset Onboard videos into the `~/VIS_Onboard` folder and the Buoy Dataset videos into the `~/Buoy` folder using the links above.
3. Set the downloaded repository folder as the MATLAB working directory. Add the video folders to the MATLAB path.

### Dependencies
- MATLAB R2018a or later.
- Required Toolboxes:
  - Image Processing Toolbox
  - Statistics and Machine Learning Toolbox
  - Computer Vision Toolbox
  - Econometrics Toolbox

## Usage
To run the code, call the function `HL_Detect_TSA` with the input video name and the number of frames, `N`, for the listener block.

#### Example:

- **For the Singapore Maritime Dataset video**  
  `vid_name = 'MVI_0788_VIS_OB.avi'`, set `N = 60`.

- **For the Buoy Dataset video**  
  `vid_name = 'buoyGT_2_5_3_4.avi'`, set `N = 20`.

### Expected Output
When you run the script, the output will be a set of annotated video frames where:
- The **parallelogram ROI** will be highlighted in yellow.
- The **detected horizon line** will be shown in red.

Additionally, the function will return a matrix containing the estimated horizon line state for each frame.

## Reproducibility
To reproduce the results from the paper:
1. Download and set up the required datasets in the corresponding folders (`~/VIS_Onboard`, `~/Buoy`).
2. Ensure the required MATLAB toolboxes are installed.
3. Run the `HL_Detect_TSA.m` script with the appropriate dataset and parameters.
4. The resulting horizon line state matrix will be saved for further analysis, and annotated frames will be displayed.

## Citation
If you use this code in your research, please cite:
@article{Agaoglu2024, title={Dynamic Region of Interest Generation for Maritime Horizon Line Detection using Time Series Analysis}, author={Ahmet Agaoglu and Nezih Topaloglu}, journal={Under Review}, year={2024}, doi={} }


## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
