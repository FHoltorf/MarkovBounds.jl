# [Background on Moment Bounding Schemes](@id background)

## Polynomial Jump-Diffusion Processes
A jump-diffusion process is dynamical system combining a deterministic evolution of the system state, called drift, with a stochastic vibrations driven by a Brownian Motion, called diffusion, and another stochastic contribution modeling discrete changes driven by Poisson counters, called jumps. The evolution of the process state ``x_t`` over time ``t`` through its state space ``X \subset \mathbb{R}^n`` is governed by the following stochastic differential equation
```math
    dx_t = f(x_t) \, dt + g(x_t) \, dW_t + \sum_{i=1}^{n_R} h(x_t) \, dN_{a_i(x_t)}
```
where $W_t$ denotes a standard ``\mathbb{R}^m``-Brownian motion and ``N_{a_i(x_t)}`` a standard Poisson counter with rate ``a_i``. The problem data is considered

* drift coefficient ``f:\mathbb{R}^n \to \mathbb{R}^n``
* diffusion matrix ``gg^\top :  \mathbb{R}^n \to \mathbb{R}^{n\times n}`` (or diffusion coefficient $ g:\mathbb{R}^n \to \mathbb{R}^{n \times m} $)
* arrival rates $a_i : \mathbb{R}^n \to \mathbb{R}$, $i  = 1,\dots, n_R$
* jumps $h_i:\mathbb{R}^n \to \mathbb{R}^n$, $i  = 1,\dots, n_R$
* state space $X$ 

A fundamental assumption in MarkovBounds.jl (and to a large extent moment bounding schemes inherently) is that the data of the jump-diffusion process under investigation can be fully characterized in terms of polynomials, i.e., all functions listed above are polynomials (component-wise) and the state space is (or can at least be outer approximated by) a [basic closed semialgebraic set](https://www.mit.edu/~parrilo/cdc03_workshop/10_positivstellensatz_2003_12_07_02_screen.pdf). Throughout, we will refer to processes that satisfies this assumption as polynomial jump-diffusion processes. A wide range of problems, in particular in the realm of stochastic chemical kinetics, lend themselves to be modeled in terms of polynomial jump-diffusion processes; if this assumption, however, is not satisfied, a simple but limited workaround is to find approximations of the data in terms of polynomials and apply the moment bounding scheme in a second step.

## Moment Bounding Schemes
The core idea behind moment bounding schemes is rather simple. But to explain it, we first need to establish some notation: Let ``y_i(t) = \mathbb{E} \left[ \prod_{k=1}^n x_k(t)^{i_k} \right]`` denote the moment corresponding to the multi-index ``i \in \mathbb{N}_{0}^n`` of a polynomial jump-diffusion process as defined the previous section. Similarly, let ``\mathbf{y}_{q}(t)``  be the truncated sequence of all multivariate moments of the process up to order ``q \in \mathbb{N}``, i.e., ``\mathbf{y}_q(t) = \{ y_i(t) | |i| \leq q \}``. Due to the notorious moment closure problem, ``\mathbf{y}_q(t)`` cannot in general be computed directly via simple simulation. To circumvent this issue, moment bounding schemes seek to identify a proxy for ``\mathbf{y}_q(t)``, say ``\tilde{\mathbf{y}}_q(t)``, which minimizes (or maximizes if upper bound is sought) a certain statistic of the process under investigation, while ensuring that ``\tilde{\mathbf{y}}_q(t)`` remains in certain ways consistent with the process under investigation (we will see shortly what that means concretely). Slightly more formally, we seek to solve an optimization problem of the form
```math
\begin{aligned} 
    \inf_{\tilde{\mathbf{y}}_q} \quad &\int_{0}^T l^\top \tilde{\mathbf{y}}_q(t) \, dt + m^\top \tilde{\mathbf{y}}_q(T) \\
    \text{s.t.} \quad & \tilde{\mathbf{y}}_q \text{ satisfies necessary consistency conditions.} 
\end{aligned}
```
The key insight underpinning all moment bounding schemes now is that a suitable choice of "necessary consistency conditions" turns the above "pseudo" optimization problem into a convex optimization problem known as generalized moment problem. The practical value of this observations lies in the fact that strong convex relaxations of these generalized moment problems are easily constructed and they can be readily solved with off-the-shelve semidefinite programming (SDP) solvers such as Mosek, SeDuMi or SDPT3. 

But what are these "necessary consistency conditions"? We won't answer this question in detail here but provide some examples and intuition for their nature. The above mentioned consistency conditions can be loosely classified as a) reflecting the dynamics of the underlying process and b) the support of its distribution. Conditions of type a) are affine relations that the moments of process have to satisfy. To derive these conditions, note that the (extended) infinitesimal generator of the process
```math
    \begin{aligned}
        \mathcal{A} : w(t,z) \mapsto &\lim_{h\to 0^+} \frac{\mathbb{E_z[w(t + h,x(t + h))]} - w(t,z)}{h} \\
                                &= \frac{ \partial w(t,z) }{\partial t } + f(z)^\top \nabla_z w(t,z) + \text{Tr}\left(gg^\top(z) \nabla_z^2 w(t,z) \right) + \sum_{i=1}^{n_R} a_i(z) w(t, h_i(z))
    \end{aligned}
```
maps polynomials to polynomials under the assumption of a polynomial jump diffusion process. The moments of a polynomial jump-diffusion process accordingly follow linear, albeit generally underdetermined, dynamics:
```math
    \frac{d}{dt}\mathbb{E}\left[\prod_{k=1}^n x_k^{i_k}(t) \right] = \mathbb{E}\left[ \mathcal{A}\prod_{k=1}^n x_k^{i_k}(t) \right] \iff \frac{dy_i}{dt}(t) = a_i^\top \mathbf{y}_q(t)
```
From this point we could in principle derive more tractable conditions; however, since this requires introduction of more technical terms and we believe it does not contribute to building better intuition, we will refer the interested reader to the references listed below instead of going through the construction here.

Conditions of type b) impose positive semidefiniteness of certain moment matrices. To understand why such conditions are in fact somehow natural, consider a one-dimensional process ``x(t)`` and a vector of the monomial basis ``b(x) = [1, x, x^2, \dots, x^d]``. Then clearly,  the moment matrix
```math
    \mathbb{E}\left[ b(x(t)) b(x(t))^\top \right] = \begin{bmatrix} 1 & y_1(t) & y_2(t)&\cdots & y_d(t) \\
															    y_1(t) & y_2(t) & y_3(t) & \cdots & y_{d+1}(t) \\
															    y_2(t) & y_3(t)  & \ddots & & y_{d+2}(t) \\
															    \vdots & \vdots & &  &\vdots \\
															    y_d(t) & y_{d+1}(t) & y_{d+2}(t) & \cdots & y_{2d}(t) 
															    \end{bmatrix}
```
must be positive semidefinite as the left-hand-side is. This argument generalizes immediately to the multivariate case. The condition,
```math
    \mathbb{E}\left[ b(x(t)) b(x(t))^\top \right]= \int_X b(x) b(x)^\top \, dP(x,t) \succeq 0,
```
can be viewed as reflecting non-negativity of the probability measure ``P(\cdot,t)`` describing the distribution of the process state at time $t$. With this intuition in mind, it follows further that for any polynomial ``p`` which is non-negative on the state space ``X``, the condition
```math
    \mathbb{E}[p(x(t))b(x(t)) b(x(t))^\top] \succeq 0
```
must hold, reflecting the support of the probability distribution ``P(\cdot,t)`` on the state space $X$. Further observe that conditions of this form translate directly into [Linear Matrix Inequalities](https://en.wikipedia.org/wiki/Linear_matrix_inequality) on the moments of the process, suggesting why the resulting problems can be tackled via SDP. 

For more details and technicalities on moment bounding schemes, please consult the references below.

## References
[1] Holtorf, Flemming, and Paul I. Barton. "Tighter bounds on transient moments of stochastic chemical systems." arXiv preprint arXiv:2104.01309 (2021).

[2] Kuntz, Juan, et al. "Bounding the stationary distributions of the chemical master equation via mathematical programming." The Journal of chemical physics 151.3 (2019): 034109.

[3] Dowdy, Garrett R., and Paul I. Barton. "Dynamic bounds on stochastic chemical kinetic systems using semidefinite programming." The Journal of chemical physics 149.7 (2018): 074103.

[4] Dowdy, Garrett R., and Paul I. Barton. "Bounds on stochastic chemical kinetic systems at steady state." The Journal of chemical physics 148.8 (2018): 084106.

[5] Sakurai, Yuta, and Yutaka Hori. "Bounding transient moments of stochastic chemical reactions." IEEE Control Systems Letters 3.2 (2018): 290-295.

[6] Sakurai, Yuta, and Yutaka Hori. "Optimization-based synthesis of stochastic biocircuits with statistical specifications." Journal of The Royal Society Interface 15.138 (2018): 20170709.

[7] Sakurai, Yuta, and Yutaka Hori. "A convex approach to steady state moment analysis for stochastic chemical reactions." 2017 IEEE 56th Annual Conference on Decision and Control (CDC). IEEE, 2017.

[8] Ghusinga, Khem Raj, et al. "Exact lower and upper bounds on stationary moments in stochastic biochemical systems." Physical biology 14.4 (2017): 04LT01

[9] Kuntz, Juan, et al. "Bounding stationary averages of polynomial diffusions via semidefinite programming." SIAM Journal on Scientific Computing 38.6 (2016): A3891-A3920.

[10] Lasserre, Jean B. Moments, positive polynomials and their applications. Vol. 1. World Scientific, 2009.

[11] Savorgnan, Carlo, Jean B. Lasserre, and Moritz Diehl. "Discrete-time stochastic optimal control via occupation measures and moment relaxations." Proceedings of the 48h IEEE Conference on Decision and Control (CDC) held jointly with 2009 28th Chinese Control Conference. IEEE, 2009.

[12] Lasserre, Jean B., Tomas Prieto‐Rumeau, and Mihail Zervos. "Pricing a class of exotic options via moments and SDP relaxations." Mathematical Finance 16.3 (2006): 469-494.

[13] Lasserre, Jean B. "Global optimization with polynomials and the problem of moments." SIAM Journal on optimization 11.3 (2001): 796-817.

[14] Helmes, Kurt, Stefan Röhl, and Richard H. Stockbridge. "Computing moments of the exit time distribution for Markov processes by linear programming." Operations Research 49.4 (2001): 516-530.

[15] Schwerer, Elizabeth. "A linear programming approach to the steady-state analysis of reflected Brownian motion." (2001): 341-368.

[16] Bhatt, Abhay G., and Vivek S. Borkar. "Occupation measures for controlled Markov processes: Characterization and optimality." The Annals of Probability (1996): 1531-1562.