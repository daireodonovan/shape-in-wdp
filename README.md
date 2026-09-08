# shape-in-wdp

A GitHub repository with the code corresponding to the paper "Shape considerations in wave-driven propulsion" by Daire O'Donovan and Graham P. Benham

## Layout of the repository

The repository is divided into subdirectories which we will discuss in alphabetical order. 

## analytical_informed

Here we compare the pure numerical optimisation to an analytically informed optimisation imposing the $\mathcal{L}_{2D}\hat{Q}=0$ condition. The reference data is stored in the data folder. 

### `param_code.jl`
Function for start-up ($v=0$) analytically informed optimisation over a square to correspond to `subcritical_square.jl`.

### `startup_compare.ipynb`
Notebook using `param_code.jl` compare start-up optimisation for Section 3

### `sub_compare.ipynb`
Notebook using `subcrit_circle.jl` to compare subcritical analytical optimisation for a circle for Section 4.1.

### `subcrit_circle.jl`
Function for general subcritical velocity ($v<1$) analytically informed optimisation for a circle.

### `sup_compare.ipynb`
Notebook using `supcrit_circle.jl` to compare supercritical analytical optimisation for a circle for Section 5.

### `supcrit_circle.jl`
Function for general supercritical velocity ($v>1$) analytically informed optimisation for a circle.
There is functionality to change the mode which was used to check with route of the analytically informed (depending on whether $\kappa$ makes real or imaginary parts in the exponential) where the two best working modes are used in `sup_compare.ipynb`.

---

## aspect_ratio

Focused on studying the effect of aspect ratio on the thrust and efficiency outputs of wave-driven propulsion (WDP).
Reference data is stored in `data` directory and relevant context for each can be found in the notebooks.

### `aspect_ratio_runfile.jl`
Julia script to fill matrix of values for different aspect ratios and velocities and saving to JLD2 file.

### `full_v_run.jl`
Julia script to fill matrix of values for 3 different aspect ratios and velocities over both sub and supercritical velocities and saving to JLD2 file.

### `full_velocity.ipynb`
Plotting data from `full_v_run.jl` output for figure 7. Note the different JLD2 files in use in an aim to increase the grid density for the narrow body.


### `sub_aspect_sweep.ipynb`
Plotting the data output from `aspect_ratio_runfile.jl` as contour plots for figure 5.


---

## compare_gen_shape

Comparing the results from including the shape of the source as a control variable against standard circular sources.

### `extend_sub.jl`
Function to take the optimisation output and plot it on a larger domain to show more wavelengths for presentation purposes. This is quicker than making a larger domain to optimise over.

### `runfile.jl`
Julia script to run each iteration of the circle and general shape optimisation for 4 different velocities and save as JLD2 files which are in the directory. 

### `subcritical_plot.ipynb`
Take circle and general shape data in the subcritical case and plot each for figures 5 and 8. The thrust and efficiency data is printed for comparison.

### `sup_plot.ipynb`
Same as `subcritical_plot.ipynb` but for the supercritical case.

## functions 
Contains frequently used functions that are used throughout the repository.

### `finite_functions.jl`
Contains the finite difference approximations of derivatives both low order and high order (high order was not used in the optimisation because it slowed it down totally).

### `h_o_optimisation.jl`
Optimisation using higher order finite difference approximation, was not used due to very slow convergence.

### `integral_functions.jl`
Functions for the integrals approximated by trapezoid rule.

## subcritical_elliptical

Directory for analysing the subcritical optimisation for an elliptical source.
Reference data is stored in the `data` directory.

### `checking_LQ.ipynb`
Notebook to run an optimisation and check $\mathcal{L}_{2D}\hat{Q}$ as well as import reference data and check prepared data. Runs a loop to check the condition over general supercritical velocities as well as doing this for reference data for different circle radius.

### `sub_analysis.ipynb`
Example notebook for plotting the result and making animations.

### `sub_elliptical.jl`
Code for subcritical elliptical optimisation.

### `validity_check.ipynb`
Notebook to check the maximum gradient scaled with maximum amplitude to test validity of the shallow water model.


## subcritical_square
Directory for analysing the subcritical optimisation for a square source. Only used for the start-up case in the paper.

### `ana_num_compare.ipynb`
Comparing the result from the numerical scheme with general subcritical boundary conditions with the Green's function solution (result plotted in figure 3).

### `numerical_solver.ipynb`
Running numerical optimisation for the start-up square source and checking $\mathcal{L}_{2D}\hat{Q}$. 

### `subcrit_solver.jl`
Subcritical square optimisation. Used for the start-up case in the paper, but can also be used for general velocities.


## supercritical_elliptical

Directory for analysing the supercritical optimisation for an elliptical source.
Reference data is stored in the `data` directory.

### `check_LQ.ipynb`
Notebook to run an optimisation and check $\mathcal{L}_{2D}\hat{Q}$ as well as import reference data and check prepared data. Runs a loop to check the condition over general subcritical velocities as well as doing this for reference data.

### `plotting.ipynb`
Example run for the supercritical optimisation and plotting. 


### `sup_elliptical.jl`
Supercritical optimisation for an elliptical source.


### `validity_check.ipynb`
Notebook to check the maximum gradient scaled with maximum amplitude to test validity of the shallow water model.

---
## `eps_contours.jl`
Julia script to output .eps plots rather than .pdf for the contour plots.

## `running_LQ.jl`
Julia script to evaluate $\mathcal{L}_{2D}\hat{Q}$ over subcritical and supercritical velocities and saving them as JLD2 files for analysis in their respective directories.


---
## To do 

- Tidy up the code so that it follows the same output convention for each function. In the development phase, I was just adding arguments when I needed them.
- Tidy up subcritical and supercritical elliptical repos since `compare_gen_shape` has the plots used in the paper.