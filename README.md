# QuantumJulia

*This README was generated automatically via DrWatson and GitHub Copilot.*

## Table of Contents
1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Reproducing the Project](#reproducing-the-project)
4. [Script Usage](#script-usage)

## Introduction

This code base is using the [Julia Language](https://julialang.org/) and
[DrWatson](https://juliadynamics.github.io/DrWatson.jl/stable/)
to make a reproducible scientific project named
> QuantumJulia

It is authored by JanSuntajs.

## Installation

To install Julia, follow the instructions in the [Julia installation manual](https://docs.julialang.org/en/v1/manual/installation/). The recommended way to install Julia is by using the Julia version multiplexer [juliaup](https://github.com/JuliaLang/juliaup).

### Getting Started

- [Getting started](https://docs.julialang.org/en/v1/manual/getting-started/)
- [Several Tutorials](https://julialang.org/learning/tutorials/)
- [Julia Tutorial for Science and Engineering](https://www.matecdev.com/posts/julia-tutorial-science-engineering.html)
- [Numerical Computing](https://www.matecdev.com/posts/julia-numerical-computing.html)
- [Zero2Hero videos (by the author of DrWatson)](https://www.youtube.com/watch?v=Fi7Pf2NveH0)

## Reproducing the Project

To (locally) reproduce this project, do the following:

0. Download this code base. You can clone the repository using git:
   ```
   git clone https://github.com/yourusername/QuantumJulia.git
   ```
   Notice that raw data are typically not included in the git-history and may need to be downloaded independently.
1. Open a Julia console and do:
   ```
   julia> using Pkg
   julia> Pkg.add("DrWatson") # install globally, for using `quickactivate`
   julia> Pkg.activate("path/to/this/project")
   julia> Pkg.instantiate()
   ```

This will install all necessary packages for you to be able to run the scripts and
everything should work out of the box, including correctly finding local paths.

## Script Usage

You may notice that most scripts start with the commands:
```julia
using DrWatson
@quickactivate "QuantumJulia"
```
which auto-activates the project and enable local path handling from DrWatson.

To run scripts in the proper environment, use:
```
julia --project=. <script_name>.jl
```

## Disclaimers

### Code Organization

Users should strive to make their code modular and avoid repetition. The `src` and `scripts` folders have distinct purposes:

- If including `file.jl` produces any output (data files, plots, or console output), it should be placed in `scripts`.
- If it contains functionality used across multiple files or pipelines, it should be in `src`.
- `src` should only contain files that define functions or types but do not produce any output.

### Project Scope

This project is designed as a scientific project using DrWatson. It is not intended to be developed as a package. For information on developing a Julia package, refer to [Developing Packages](https://julialang.org/contribute/developing_package/).
```
