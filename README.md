# Spectral Computations

This repository accompanies the book *Infinite-Dimensional Spectral Computations* by Matthew J. Colbrook. It contains code examples, supporting routines, exercise-solution code, Koopman-algorithm examples, and data used throughout the chapters.

## Repository Structure

The repository is organised as follows:

- `chapter_1/`, `chapter_3/`, ..., `chapter_11/` — Code and examples corresponding to chapters of the book
- `utils/` — Shared functions and utilities used across multiple chapters
- `Koopman_algorithms/` — Code related to Koopman algorithms and associated examples
- `exercise_sol_code/` — Code for solutions to selected exercises from the book
- `data_online/` — Smaller data files included directly in the repository

The chapter folders are largely self-contained, but many examples rely on functions provided in the `utils/` directory. Please add the relevant folders to your MATLAB path before running the examples.

## Getting Started

To run the code, first clone the repository:

    git clone https://github.com/MColbrook/spectralcomputations.git
    cd spectralcomputations

Then add the relevant folders to your MATLAB path. In particular, many examples require `utils/` to be on the path. Some examples may also require chapter-specific folders, `Koopman_algorithms/`, `exercise_sol_code/`, or the relevant data directories.

## External Dependencies

Some of the code in this repository uses Chebfun.

You can download Chebfun from:

[Chebfun](https://www.chebfun.org/)

Please ensure that Chebfun is installed and added to your MATLAB path before running examples that depend on it.

## Data Files

Some examples require larger data files that are not stored directly in this repository.

These can be downloaded from:
[Book_data_sets](https://www.dropbox.com/scl/fo/w8j7pmegivi37ngi6c4u9/AO1rgz5HBct8_3-qnEtA6EE?rlkey=zjcycyn2u7s02pin6lrh6usil&st=45v8z4d3&dl=0)

After downloading the data files:

- Place them in a suitable directory on your system
- Add this directory to your MATLAB path for the relevant examples
- Alternatively, update the file paths in the relevant scripts so that they point to the downloaded data

The folder `data_online/` contains data files that are small enough to be included directly in this repository.

## License

This repository is released under the MIT License. See `LICENSE` for details.

## Contact

For questions or issues, please open an issue on GitHub or visit:

[Matthew J. Colbrook's webpage](https://www.damtp.cam.ac.uk/user/mjc249/home.html)
