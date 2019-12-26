# Gradient-based boundary (GBB)

- general description

## General idea
## Prerequisites
## Comments


## necessary packages
os
random
numpy
nibabel
matplotlib
one function from lib
scipy
sys
cv2

Algorithm
- step 1: choose randomly vertex point on surface mesh
- step 2: get gradient along one direction
- step 3: update neighborhood
- step 4: update cost function (BBR)
- step 5: stop if the maximum number of iterations has been reached or if cost function reached steady-state; otherwise modify $L(t)$ and $\Theta(t)$ and continue with step 1

### Update neighborhood
- $x(t+1)=x(t)+\Theta(t)\cdot L(t)(y(t)-x(t))$
- $t$: time step
- $y(t)$: proposed new position
- $\Theta(t)$: neighborhood kernel (coupling strength between neighbors)
- $L(t)$: learning rate
- $\lambda$: time constant
- $\sigma$: neighborhood size
- $\Theta(t)_{ci}=\exp(-\|r_c-r_i\|^2/2\sigma^2(t))$; $\sigma(t)=\sigma_0\exp(-t/\lambda)$
- $L(t)=L_0\exp(-t/\lambda)$
- $L(t)=A/(t+B)$ (option 2)
- decrease $\sigma$ and $L(t)$ over time

### Comments
- a “rule” of thumb is that, for good statistical accuracy, the number of steps must be at least 500 times the number of network units
- for approximately the first 1000 steps, $L(t)$ should start with a value that is close to unity, thereafter decreasing monotonically
- we can start with a smaller learning rate from the beginning
- decrease neighborhood size over time
- at the end: $\Theta(t)=0$ (important for convergence) and $L(t)=\text{small}$
- decrease is essential that also steps further back still have influence on the refinement
- neighborhood function should not decay too fast
- start with large neighborhood kernel: e.g. 10 mm fwhm
