# MendelLocationScores

This [Julia](http://julialang.org/) package maps a trait via the method of Location Scores, i.e., multipoint linkage analysis. MendelLocationScores is one component of the umbrella [OpenMendel](https://openmendel.github.io) project.

[![](https://img.shields.io/badge/docs-current-blue.svg)](https://OpenMendel.github.io/MendelLocationScores.jl)

## Installation

*Note: The three OpenMendel packages (1) [SnpArrays](https://openmendel.github.io/SnpArrays.jl/latest/), (2) [Search](https://openmendel.github.io/Search.jl), and (3) [MendelBase](https://openmendel.github.io/MendelBase.jl) must be installed before any other OpenMendel package will run. It is easiest if these three packages are installed in the above order and before any other OpenMendel package.*

Within Julia, use the package manager to install MendelLocationScores:

    Pkg.clone("https://github.com/OpenMendel/MendelLocationScores.jl.git")

This package supports Julia v0.4 and v0.5.

## Data Files

To run this analysis package you will need to prepare a Control file and have your data files available. The Control file holds the names of your data files and any optional parameters for the analysis. Details on the general format and contents of the Control and data files can be found on the MendelBase [documentation page](https://openmendel.github.io/MendelBase.jl). Descriptions of the specific options available within the MendelLocationScores analysis package are in its [documentation page](https://openmendel.github.io/MendelLocationScores.jl).

There are example data files in the "docs" subfolder of each Mendel package, for example, ~/.julia/v0.5/MendelLocationScores/docs.

## Running the Analysis

To run this analysis package, first launch Julia. Then load the package with the command:

     julia> using MendelLocationScores

Next, if necessary, change to the directory containing your files, for example,

     julia> cd("~/path/to/data/files/")

Finally, to run the analysis using the parameters in the control file Control_file.txt use the command:

     julia> LocationScores("Control_file.txt")

*Note: The package is called* MendelLocationScores *but the analysis function is called simply* LocationScores.

## Citation

If you use this analysis package in your research, please cite the following reference in the resulting publications:

*Lange K, Papp JC, Sinsheimer JS, Sripracha R, Zhou H, Sobel EM (2013) Mendel: The Swiss army knife of genetic analysis programs. Bioinformatics 29:1568-1570.*

<!--- ## Contributing
We welcome contributions to this Open Source project. To contribute, follow this procedure ... --->

## Acknowledgments

This project is supported by the National Institutes of Health under NIGMS awards R01GM053275 and R25GM103774 and NHGRI award R01HG006139.
