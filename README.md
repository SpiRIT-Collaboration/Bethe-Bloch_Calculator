### configuration files
- style.h: configuring drawing style
- functions.h: containing functions for global usage

### making data
1. make_pid.C         -> pid histogram
2. proton_fit.C       -> fit data (function parameter c1, c2)
3. calculate_mass.C   -> mass data tree

### Modified Bethe-Bloch Formula
> by Leo. Techniques for Nuclear and Particle Physics Experiments

<img src="figures/formula1.png" alt="drawing" width="400px"/>
with
<img src="figures/formula2.png" alt="drawing" width="400px"/>
where Z(′) and M(′) are charge and weight of the transporting-particle(absorbing-material), and me is electron mass. The fitting parameters C1 and C2 are shared by all particle species.
