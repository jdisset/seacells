#Evolved Development Strategies of Artificial Multicellular Organisms
This repository contains the code and some of the results of an artificial life experiment in which cells must survive as long as possible in a demanding environment.
Light is present in the water, nutrients are in the ground, cells need to keep both their levels of light and nutrients above zero, they can absorb light and nutrients when exposed to them and can share these two forms of energy. 

Organisms are evolved using a standard GA (https://github.com/jdisset/gaga). Simulations start with one single cell, its dna is an artificial gene regulatory network (https://github.com/jdisset/grgen). World and cells physics are simulated using a custom cellular physics model (https://github.com/jdisset/MecaCell).

# Video captures
In the following video captures, the nutrients are represented as small spheres. Their colors indicate the available amount (green = full, red = empty), and the cells colors represent their level in light and nutrients.

Click on the images to access the videos.

## Simple novelty (Nm0)
[![Surv Only](https://raw.githubusercontent.com/jdisset/seacells/master/images/novsurv5.jpg)](https://vimeo.com/157664099)

[Same organism without ground view](https://vimeo.com/157666837)
## Multi novelty (Nm1)
[![Surv Only](https://raw.githubusercontent.com/jdisset/seacells/master/images/multinov.jpg)](https://vimeo.com/157664106)
## Capture novelty (Nm2)
[![Surv Only](https://raw.githubusercontent.com/jdisset/seacells/master/images/novcapture.jpg)](https://vimeo.com/157666811)
## No novelty, survival only (local optima toward which 70% of the runs converged )
[![Surv Only](https://raw.githubusercontent.com/jdisset/seacells/master/images/survonly.jpg)](https://vimeo.com/157664875)

This repository will be regularly updated with new videos of organisms evolved in different world configurations.


