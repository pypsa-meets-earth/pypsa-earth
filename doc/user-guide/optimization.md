# Optimization Models in Energy Systems

The aim of this section is to provide an overview of the most common optimization approaches used in energy system models.
In particular, it covers capacity expansion models and economic dispatch models, and it distinguishes between deterministic and stochastic formulations.
Before addressing these categories, we introduce some general definitions.
In general, an optimization problem can be expressed as:


   \begin{aligned}
       \min_{\mathbf{x}} \quad & f(\mathbf{x}) \quad \text{(objective function)}
   \end{aligned}

subject to:



   \begin{aligned}
       & g_i(\mathbf{x}) \le 0 \quad \forall i \quad \text{(inequality constraints)} \\
       & h_j(\mathbf{x}) = 0 \quad \forall j \quad \text{(equality constraints)}
   \end{aligned}


where:

-  $f(\mathbf{x})$ is the objective function, i.e., the function
   to be minimized (or maximized).

-  $\mathbf{x}$ are the *decision variables*, i.e., the variables
   you can choose to minimize (maximize) the objective function.

-  $g_i$ and $h_j$ represent the set of inequality and
   equality constraints.


# Example: Simple linear optimization problem

Consider the following problem:



   \begin{aligned}
       \min_{x_1, x_2} \quad & x_1 + x_2 \label{eq:obj_simple}
   \end{aligned}

subject to:



   \begin{aligned}
       & x_1 \geq 0 \quad \label{eq:const1} \\
       & x_2 \geq 0 \quad \label{eq:const2} \\
       & x_1 + 2x_2 = 4 \quad \label{eq:const3}
   \end{aligned}

**Explanation:**

-  The *objective function* aims to
   minimize the sum $x_1 + x_2$.

-  The *decision variables* are $x_1$ and $x_2$.

-  Constraints impose that both variables must be
   non-negative.

-  Constraint $x_1 + 2x_2 = 4$ requires that the combination
   of $x_1$ and $x_2$ must satisfy a given requirement
   (e.g., total production must meet demand).

**Optimal solution:** Since $x_1$ is less ‘expensive’ in the
constraint `x_1 + 2x_2 = 4` (it counts 1 unit toward the
total while $x_2$ counts 2), the optimal choice is to use as much
$x_2$ as needed to meet the constraint with minimal total
$x_1 + x_2$. Solving:



gives an objective value $0 + 2 = 2$.


# Capacity expansion problem

Let us assume a single-year timeframe, with the aim of minimizing the
*total cost* of the system under technical and policy constraints.

We consider a system composed of:

-  A photovoltaic (PV) plant with capacity $P_{\text{PV}}$ [MW].

-  A gas generator with capacity $P_{\text{gas}}$ [MW].

-  A fixed hourly electricity demand $D(t)$ [MW] for
   $t=1,\dots,T$.

The optimization problem can be formulated as:



   \begin{aligned}
       \min_{\substack{P_{\text{PV}},\, P_{\text{gas}}, \\ p_{\text{PV}}(t),\, p_{\text{gas}}(t)}}
       & C_{\text{PV}}^{\text{cap}} \cdot P_{\text{PV}} + C_{\text{gas}}^{\text{cap}} \cdot P_{\text{gas}} \nonumber \\
       & + \sum_{t=1}^T c_{\text{gas}} \cdot p_{\text{gas}}(t) \cdot \Delta t \label{eq:obj_energy}\end{aligned}

subject to:



   \begin{aligned}
       & p_{\text{PV}}(t) + p_{\text{gas}}(t) = D(t) \quad \forall t \label{eq:bal} \\
       & 0 \leq p_{\text{PV}}(t) \leq \min\{G(t) \cdot \eta_{\text{PV}},\; P_{\text{PV}}\} \quad \forall t \label{eq:pv_lim} \\
       & 0 \leq p_{\text{gas}}(t) \leq P_{\text{gas}} \quad \forall t \label{eq:gas_lim} \\
       & P_{\text{PV}} \geq 0, \quad P_{\text{gas}} \geq 0 \label{eq:cap_nonneg}\end{aligned}

**Interpretation:**

-  The first line of the objective function represents the
   *design cost* (capital expenditure), depending on installed
   capacities.

-  The second line of the objective function represents
   the *operational cost* (fuel cost for the gas generator).

-  The equality constraint enforces the energy balance, i.e. the sum of generation must
   match the demand.

- The inequality constraints enforce the technical limits. The gas production is bounded by the
  gas generator size, while the PV production is bounded by the PV size and irradiance availability.
  Both productions must be positive (generators cannot act like loads).

Since capacities are part of the decision variables, this is a capacity expansion model:
it optimizes both design and operation. Here, design refers to the sizing phase, i.e.
the choice of the optimal installed capacities for the system components,
while operation determines their temporal dispatch. A realistic example could be the evolution
of the European power system in 2050 under specific technical, economic, or environmental assumptions.


# Economic dispatch problem

If capacities are given and fixed, the problem becomes a special case of
capacity expansion known as *economic dispatch* (ED).

**Parameters**

-  $P_{\text{PV}}$ [MW]: installed PV capacity (fixed).

-  $P_{\text{gas}}$ [MW]: installed gas turbine capacity (fixed).

-  $D(t)$ [MW], $t=1,\dots,T$: electricity demand.

-  $G(t)$ [MWh], $\eta_{\text{PV}}$
    : PV availability and conversion factor.

-  $c_{\text{gas}}$ [€ / MWh], $\Delta t$ [h]: gas
   marginal cost and time step.

**Decision variables (operation)**



   p_{\text{PV}}(t) [MW], \qquad
   p_{\text{gas}}(t) [MW] \qquad \forall t=1,\dots,T

**Optimization problem**



   \begin{aligned}
       \min_{\{p_{\text{PV}}(t),\, p_{\text{gas}}(t)\}} \quad
       & \sum_{t=1}^T c_{\text{gas}} \, p_{\text{gas}}(t)\, \Delta t
       \label{eq:ed_obj}\end{aligned}

subject to



   \begin{aligned}
       & p_{\text{PV}}(t) + p_{\text{gas}}(t) = D(t)
         \quad \forall t
         && \text{(demand balance)}
         \label{eq:ed_balance} \\
       & 0 \le p_{\text{PV}}(t) \le \min\!\left\{ G(t)\,\eta_{\text{PV}},\; P_{\text{PV}} \right\}
         \quad \forall t
         && \text{(PV availability and capacity)}
         \label{eq:ed_pv_bounds} \\
       & 0 \le p_{\text{gas}}(t) \le P_{\text{gas}}
         \quad \forall t
         && \text{(gas capacity limit)}
         \label{eq:ed_gas_bounds}\end{aligned}

A realistic example could be the validation of a given model with
historical data, where capacities are set (historical ones) and only
operation is optimized.


# Stochastic Optimization

So far, we have assumed *deterministic* optimization: all input time series
(demand, solar irradiance) and parameters (natural gas price) are
perfectly known. However, in real life we often face uncertainty.

## Birthday party example.

Suppose tomorrow is your birthday and you are going to have a party. You
invited $Y=20$ people and everyone wants a pizza, but nobody has
confirmed their presence yet. This means you do not know the actual
number of guests $y$ who will show up. You must decide *today* how
many pizzas $x$ to order, at a cost of 10 € each. If more guests
arrive than pizzas ordered ($y > x$), you will need to buy extra
pizzas *last-minute* at 16 € each. There is no refund for leftovers.

If $y \le x$ (over-ordering), you spend $10x$, while
$10y$ would have been enough. If $y > x$ (under-ordering),
you spend $10x$ plus $16(y-x)$ for the extra pizzas. The
challenge: $x$ must be chosen **today** — before knowing
$y$.

Assume three equally likely scenarios:

- 0 guest, nobody is coming
- 5 guests, only your best friends are coming
- 20 guests, everyone is coming

Therefore, $y \in \{0, 5, 20\}$ guests. These scenarios have the same
probability of occurring.


**Scenario-by-scenario explanation and expected cost**

- **Case A — \(0 <= x <= 5\)** - under-ordering can already happen at y=5:

  - Scenario y=0: no guests → you only pay the pre-ordered pizzas
    \(C(x,0)=10x\).

  - Scenario y=5: if x<5, you are short of \(5-x\) pizzas → last-minute at 16€ each
    \(C(x,5)=10x+16(5-x)\). If x=5, the extra term is 0.

  - Scenario y=20: you are short of \(20-x\) pizzas → last-minute at 16€ each
    \(C(x,20)=10x+16(20-x)\).


  **Expected cost**
  Each scenario has the same probability, so its contribution is weighted 1/3:

  .. math``
 \mathbb{E}[C(x)]
 = \tfrac13\big(10x\big)
 + \tfrac13\big(10x+16(5-x)\big)
 + \tfrac13\big(10x+16(20-x)\big)
 = \frac{400 - 2x}{3}
 = 133.33 - \tfrac{2}{3}x.

``

- **Case B — \(5 < x <= 20\)** - no shortage at y=5, only at y=20:

  - Scenario y=0: no guests → you only pay the pre-ordered pizzas
    \(C(x,0)=10x\).

  - Scenario y=5: enough pizzas \(x>5\) → no last-minute purchase
    \(C(x,5)=10x\).

  - Scenario y=20: you are short of \(20-x\) pizzas → last-minute at 16€ each
    \(C(x,20)=10x+16(20-x)=320-6x\).

  **Expected cost**

  .. math``
 \mathbb{E}[C(x)]
 = \tfrac13\big(10x\big)
 + \tfrac13\big(10x\big)
 + \tfrac13\big(10x+16(20-x)\big)
 = \frac{320 + 14x}{3}
 = 106.67 + \tfrac{14}{3}x.

``

**Minimizer.** Since $\mathbb{E}[C(x)]$ is decreasing on [0,5] and increasing on [5,20],
the minimum is attained at the boundary $x^\star=5$, with



   \mathbb{E}[C(5)] = \frac{390}{3} = 130.


**Interpretation:** Ordering 5 pizzas perfectly covers the medium
scenario, avoids over-ordering in the low scenario, and limits the
expensive last-minute purchases in the high scenario. Here, $x$ is
a *here-and-now* decision taken under uncertainty, while the number of
extra pizzas (if needed) is a *wait-and-see* decision made after the
actual scenario is revealed.


## Two-stage stochastic formulation.

Let $\omega \in \Omega$ be a scenario with probability
$p_\omega$. First-stage variables $x$ are chosen “here and
now”, before knowing which scenario will occur. Second-stage variables
$y_\omega$ are chosen “wait-and-see”, after the specific scenario
$\omega$ is revealed.

The two-stage stochastic problem can be written as:



   \begin{aligned}
       \min_{x,\,\{y_\omega\}_{\omega\in\Omega}} \quad &
       f(x) + \sum_{\omega \in \Omega} p_\omega \, g(x, y_\omega, \omega) \label{eq:stoc_compact} \\
       \text{s.t.} \quad & A_\omega x + B_\omega y_\omega \ge b_\omega, \quad \forall \omega \in \Omega, \\
                         & x \ge 0, \quad y_\omega \ge 0 \quad \forall \omega \in \Omega.\end{aligned}

Here:

-  $f(x)$ represents the first-stage cost or contribution to the
   objective.

-  $g(x, y_\omega, \omega)$ represents the second-stage cost or
   contribution, depending on scenario $\omega$.

-  The constraints must hold for every scenario $\omega$.

In the birthday party analogy:

-  $x$ = pizzas ordered today (*first stage*);

-  $y_\omega$ = guests in scenario $\omega$ (*second
   stage*).

The same structure applies to energy systems: first-stage = investment
decisions (capacities), second-stage = operational decisions (dispatch)
under different scenarios of demand, renewable generation, or fuel
prices.


## Capacity expansion under uncertainty (three gas price scenarios).

We consider the same system as in the deterministic case:

-  A photovoltaic (PV) plant with capacity $P_{\text{PV}}$ [MW].

-  A gas generator with capacity $P_{\text{gas}}$ [MW].

-  A fixed hourly electricity demand $D(t)$ [MW],
   $t=1,\dots,T$.

Now, however, the natural gas price is uncertain. We define three
scenarios $\omega \in \{1,2,3\}$ with probabilities
$p_\omega$:



representing, for example, *low*, *medium*, and *high* price conditions.

**Decision structure:**

-  **First-stage (here-and-now)**: capacities $P_{\text{PV}}$,
   $P_{\text{gas}}$ (same for all scenarios).

-  **Second-stage (wait-and-see)**: dispatch profiles
   $p_{\text{PV}}(t,\omega)$, $p_{\text{gas}}(t,\omega)$
   (adapted to each scenario’s fuel cost).

**Two-stage stochastic formulation:**



   \begin{aligned}
       \min_{P_{\text{PV}},\, P_{\text{gas}}} \quad
       & C_{\text{PV}}^{\text{cap}} \cdot P_{\text{PV}} + C_{\text{gas}}^{\text{cap}} \cdot P_{\text{gas}} \nonumber \\
       & + \sum_{\omega=1}^3 p_\omega \; Q(P_{\text{PV}},P_{\text{gas}},\omega) \label{eq:stoc_obj_energy}\end{aligned}

where the second-stage operational cost for scenario $\omega$ is:



   \begin{aligned}
       Q(P_{\text{PV}},P_{\text{gas}},\omega) \;=\;
       & \sum_{t=1}^T c_{\text{gas}}^{(\omega)} \cdot p_{\text{gas}}(t,\omega) \cdot \Delta t\end{aligned}

subject to, for each scenario $\omega$:



   \begin{aligned}
       & p_{\text{PV}}(t,\omega) + p_{\text{gas}}(t,\omega) = D(t)
         && \forall t \quad \text{(demand balance)} \label{eq:stoc_bal} \\
       & 0 \le p_{\text{PV}}(t,\omega) \le \min\{ G(t) \cdot \eta_{\text{PV}},\; P_{\text{PV}}\}
         && \forall t \quad \text{(PV limit)} \label{eq:stoc_pv_lim} \\
       & 0 \le p_{\text{gas}}(t,\omega) \le P_{\text{gas}}
         && \forall t \quad \text{(gas capacity)} \label{eq:stoc_gas_lim}\end{aligned}

**Interpretation:**

-  In the deterministic version, $c_{\text{gas}}$ is known, so the
   model optimizes capacities and dispatch for that single case.

-  In the stochastic version, the model chooses capacities that minimize
   investment cost plus the *expected* operational cost across all gas
   price scenarios.

-  The optimal solution is a compromise: it may not be optimal for any
   single scenario, but it is best *on average*, given the probabilities
   $p_\omega$.

**Note:** If multiple investment stages are allowed (e.g., building some
capacity now and more in 10 years), the formulation extends to a
*multi-stage* stochastic problem, where each stage has its own set of
here-and-now decisions, followed by scenario-dependent operational
decisions.


# Sensitivity Analysis

**Sensitivity analysis** is performed *after* solving the optimization
problem, by varying one or more parameters to see how the optimal
solution changes.



   Solve a deterministic CE model, then vary the fuel cost
   $c_{\text{gas}}$ from 50 to 100 €/MWh to see how the optimal
   capacities change.



   - **Sensitivity analysis**: changes are explored *post-optimization*.
   - **Stochastic optimization**: uncertainty is included *within* the
     optimization process.



# Representative days.

When optimizing over an entire year with hourly resolution
($T=8760$), the problem size can become computationally heavy. A
common approach is to select a reduced set of *representative days* that
capture the main variability of demand and renewable generation
profiles. Each representative day $d \in \mathcal{D}$ is assigned
a *weight* $w_d$ indicating how many real days it represents in
the year.

This method is *deterministic*: there is only one scenario, but the time
domain is reduced. The weights ensure that the reduced problem still
approximates the annual cost and performance.

The optimization problem becomes:



   \begin{aligned}
       \min_{x,\,\{y_d\}_{d\in\mathcal{D}}} \quad &
       f(x) + \sum_{d \in \mathcal{D}} w_d \, g(x, y_d, d) \label{eq:repdays_obj} \\
       \text{s.t.} \quad & A_d x + B_d y_d \ge b_d, \quad \forall d \in \mathcal{D}, \\
                         & x \ge 0, \quad y_d \ge 0 \quad \forall d \in \mathcal{D}.\end{aligned}

The structure is similar to the stochastic formulation, but:

-  The $w_d$ are *weights*, not probabilities: they sum to the
   number of days in the year (e.g., 365), not to 1.

-  The problem is deterministic: all representative days are solved
   jointly as part of a single time series approximation.

-  There is no uncertainty: $d$ indexes clusters of days, not
   future scenarios.


# References

-  [Theory of convex optimization <sec-simple-linear-optimization-model>] – book with everything about convex optimization

   | S. P. Boyd, L. Vandenberghe, Convex Optimization, version 29 Edition, Cambridge University Press.

-  [Theory of stochastic programming <sec-stochastic-optimization>] - book with focus on stochastic programming

   | G. Infanger, Planning under uncertainty: solving large-scale stochastic linear programs, UNT Digital Library, California, 1992, report accessed on December 9, 2024. URL https://digital.library.unt.edu/ark:/67531/metadc1114558/

-  [Application of stochastic programming to design/operation optimization problems + theory explanation <sec-two-stage-formulation>]

   | G. Mavromatidis, K. Orehounig, J. Carmeliet, Design of distributed energy systems under uncertainty: A two-stage stochastic programming approach, Applied Energy 222 (2018) 932–950.   doi:https://doi.org/10.1016/j.apenergy.2018.04.019. URL https://www.sciencedirect.com/science/article/pii/S0306261918305580

   | H. Teichgraeber, A. R. Brandt, Optimal design of an electricity-intensive industrial facility subject to electricity price uncertainty: Stochastic optimization and scenario reduction, Chemical Engineering Research and Design 163 (2020) 204–216. doi:https://doi.org/10.1016/j.cherd.2020.08.022. URL https://www.sciencedirect.com/science/article/pii/S026387622030441X

-  [Game theory for energy systems <sec-two-stage-formulation>] - course explaining stochastic programming in energy systems, with focus on electricity markets

   | J. Kazempour, Advanced optimization and game theory for energy systems - youtube. URL https://www.youtube.com/playlist?list=PLe7H9pun_r8YHoGv0TnYxUsgbj0xAJmMR

-  [Applications on PyPSA <sec-capacity-expansion-problem>] - capacity expansion and economic dispatch models

   | C. Gallego-Castillo, M. Victoria, PyPSA-Spain: An extension of PyPSA-Eur to model the Spanish energy system 60 101764. doi:10.1016/j.esr.2025.101764. URL https://www.sciencedirect.com/science/article/pii/S2211467X25001270

   | K. Kwak, W. Son, Y. Yang, J. Woo, PyPSA-Korea: An open-source energy system model for planning Korea’s sustainable energy transition 13 5677–5691. doi:10.1016/j.egyr.2025.05.018. URL https://www.sciencedirect.com/science/article/pii/S2352484725002963
