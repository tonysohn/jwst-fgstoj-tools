# jwst-fgstoj-tools
This repository serves as a central place for maintaining tools related to the JWST FGS-toJ frame calibration.

- **`FGStoJ_Calibration.py`**
This tool was directly translated from the matlab tool `FGStoJ_Calibration.m` provided by NGAS. This tool does NOT solve the full matrix equation. Instead, it fixes the clocking angle for FGS1 and calculates the angular difference between the (x_FGSidl, y_FGSidl) coordinates predicted using the old FGS-to-J matrix and the observed positions. This angular difference is applied to the old matrix to come up with a new matrix. For this reason, the tool works with any number of calibration sources.

- **`FGStoJ_Calibration_FullSolution.py`**
This tool is for solving the matrix equation in its full form - see the accompanying PDF file for details on how the matrix equation is set up and how to solve it. The solution requires at least 3 calibration sources.

- **`vtoj_siaf.py`**
This tool was written to calculate the `V2Ref, V3Ref, V3IdlYAngle` SIAF parameters for the psudo-aperture J-FRAME for FGS1. Everytime the FGS-to-J matrix is updated, this tool must be run to update these SIAF parameters for the J-FRAME aperture. The `vtoj_siaf.ipynb` Jupyter notebook provides examples on how to run the tool. 

- **`vtoj_siaf_rotation.py`**
This is an alternative tool to the `vtoj_siaf.py` using the `Rotation` method in `scipy.spatial.transform` package. The tool was written to verify `vtoj_siaf.py` and results in slightly different outputs. DO NOT USE this tool for final results that go into the SIAF updates - use `vtoj_siaf.py` instead.
