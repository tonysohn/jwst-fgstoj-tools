# jwst-fgstoj-tools
This repository serves as a central place for maintaining tools related to the JWST FGS-toJ frame calibration. Below are the high-level descriptions on what each tool does and how to use them.

- **`fgstoj_calibration.py`**
This code includes tools for solving the matrix equation to obtain the FGSics to J-frame matrix. 

    (1) The method provided by NGAS in their matlab tool `FGStoJ_Calibration.m`. This tool does NOT solve the full matrix equation. Instead, it calculates the angular difference between the (x_FGSidl, y_FGSidl) coordinates predicted using the old FGS-to-J matrix and the observed positions. This angular difference is applied to the old matrix to come up with a new matrix. For this reason, the tool works with any number of calibration sources.
    
    (2) This tool is for solving the matrix equation in its full form - see the accompanying PDF file for details on how the matrix equation is set up and how to solve it. The solution requires at least 3 calibration sources.


- **`vtoj_siaf.py`**
This tool was written to calculate the `V2Ref, V3Ref, V3IdlYAngle` SIAF parameters for the psudo-aperture J-FRAME for FGS1. Everytime the FGS-to-J matrix is updated, this tool must be run to update these SIAF parameters for the J-FRAME aperture. The `vtoj_siaf.ipynb` Jupyter notebook provides examples on how to run the tool.


- **`fgs2tofgs1_calibration.py`**
This tool is for calculating the FGS2 to FGS1 matrix. Required inputs are SIAF alignment parameters \[V2Ref, V3Ref, V3IdlYAngle\] (must be in this order) for FGS2 and FGS1.
