# Impulses and Physiological States in Theoretical Models of Nerve Membrane
## Richard Fitzhugh

This article presented in $1961$ by Richard Fitzhugh [1] proposes a revision and further study of the Van der Pol's equation for a relaxation oscillator by the addition of terms to produce a pair of non-linear differential equations with either a stable singular point or a limit cycle. Unlike the Hodgkin - Huxley model (HH model), the resulting "BVP model" has only two variables of state, representing excitability and refractoriness. Both models are studied and framed as members of a wider class of excitable-oscillatory systems. An in-depth analysis on the phase plane is carried out with the aim to highlight physiological states of a nerve fiber (resting, active, refractory, enhanced, depressed, etc.) to form a "physiological state diagram". 
Furthermore, a properly chosen projection from the $4$-dimensional HH phase space [2] onto a plane produces a very similar diagram, showing the underlying relationship between the two models. Impulse trains occur in the BVP and HH models for a range of constant applied currents which make the singular point representing the resting state unstable.
I studied and reproduced many of the properties examined by R. Fitzhugh and found a way to produce a finite train of impulses. Then, I highlighted an Hopf transition [3] for the case of the unstable point of the limit cycle ($z=-0.4$). In the end, I implemented the Gillespie algorithm to simulate the associated microscopic stochastic model.

All the code with numerical integration, plots and the Gillespie algorithm can be found for further reference at: [GitHub/RiccardoTancredi/Project](https://github.com/RiccardoTancredi/P.M.L.S./tree/main/Project).

Slides at: [RiccardoTancredi/GoogleSlides](https://docs.google.com/presentation/d/1_MTUj5nJrZAfwSivjQhx4RRLk9Se5WQFmq35NNqR35E/edit?usp=sharing)


### References
[1] R. FitzHugh, Impulses and physiological states in theoretical models of nerve membrane, Biophysical journal, 1, 445 (1961).

[2] R. Fitzhugh, Thresholds and plateaus in the hodgkin-huxley nerve equations, The Journal of general physiology 43, 867 (1960).

[3] E. Wallace, M. Benayoun, W. Van Drongelen, and J. D. Cowan, Emergent oscillations in networks of stochastic spiking neurons, Plos one 6, e14804 (2011).


