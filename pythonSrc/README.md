# Python translation of Fortran Shuttle Project
# README for pythonSrc

This folder contains the Python translation of the Fortran Shuttle Project. The main script is `shuttle_filter.py`.

## How to Run

1. Ensure you have Python 3 installed.
2. Make sure the input data file `FILTER_IN.DAT` is present in the `fortranSrc` folder (or update the path in the script).
3. Run the script:

```bash
python3 pythonSrc/shuttle_filter.py
```

## Features
- Console-based NASA shuttle ASCII art
- Reads filter data from the Fortran input file
- Prompts for gamma value
- Outputs (W, Iterations) table
- Displays an ASCII scatter plot of Iterations vs W (no external plotting libraries required)

---
Original Fortran by Steve Winward. Python conversion by Copilot.
