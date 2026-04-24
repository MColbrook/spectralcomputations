# Spectral Computations

This repository accompanies the book *Infinite-Dimensional Spectral Computations* by Matthew J. Colbrook. It contains code examples, supporting routines, and data used throughout the chapters.

## Repository Structure

The repository is organised by chapter:

- `chapter1/`, `chapter2/`, … — Code and examples corresponding to each chapter of the book  
- `routines/` — Shared functions and utilities used across multiple chapters  

The chapter folders are largely self-contained, but many examples rely on functions provided in the `routines/` directory.

## Getting Started

To run the code:

1. Clone the repository:

    git clone https://github.com/MColbrook/spectralcomputations.git  
    cd spectralcomputations

2. Add the relevant folders to your path (MATLAB or Python, depending on your setup). In particular, ensure that the `routines/` directory is on your path, as many examples depend on it.

## External Dependencies

Some of the code in this repository uses Chebfun.

You can download Chebfun from:  
https://www.chebfun.org/

Please ensure Chebfun is installed and added to your MATLAB path before running examples that depend on it.

## Data Files

Some examples require larger data files that are not stored directly in this repository.

These can be downloaded from:  
[Book_data_sets](https://www.dropbox.com/scl/fo/w8j7pmegivi37ngi6c4u9/AO1rgz5HBct8_3-qnEtA6EE?rlkey=zjcycyn2u7s02pin6lrh6usil&st=45v8z4d3&dl=0)

After downloading:
- Place the data files in a suitable directory on your system  
- Add this directory to your path so that the relevant examples can locate the data  

## Contact

For questions or issues, please open an issue on GitHub or visit https://www.damtp.cam.ac.uk/user/mjc249/home.html
