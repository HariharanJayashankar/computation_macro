---
title: Menu Cost Model
author: Hariharan Jayashankar
---

# Timing

Within a period $t$

1. Idiosyncratic and aggregate shocks realize
2. Firms make pricing decisions
3. Aggregate prices realize

# Final Good Aggregators

There's a final goods producer who is just a CES aggregator. They simply exist to give us three things:
- The demand for each intermediate good
- The Aggregate price index
- Aggregate output

Aggregate output is

$$
Y_t = \left( \int y_{it}^{\frac{\epsilon - 1}{\epsilon}} \right)^{\frac{\epsilon}{\epsilon - 1}}di
$$

Demand for a variety $i$ is given by

$$
y_{it} = p_{it}^{-\epsilon} Y_t
$$

Where $p_{it}$ is the relative price of variety $i$. (i.e. their nominal price divided by the aggregate price index)

The aggregate price index written in terms of relative prices implies

$$
\left(\int p_{it}^{1 - \epsilon} di \right)^{\frac{1}{1-\epsilon}} = 1
$$


# Intermediate Goods

The production function is given by $y_{it} = e^{Z_t a_{it}} L_{it}$

Real operating profit is given by:

$$
\Pi_{it} = p_{it} y_{it} - w_t L_t
$$

Where $w_t$ is the real wage.

Substituting in the demand for good $i$ from the final goods aggregor we get


$$
\Pi(p_{it}, a_{it}, Z_t) = \left(p_{it}^{1 - \epsilon} - p_{it}^{-\epsilon} \frac{w_t}{e^{Z_t a_{it}}} \right) Y_t
$$

The recursive representation of the intermediate goods producer is given by:

$$
V^{NA}_t(p_{i}, a_{i}) = \Pi(\frac{p_{i}}{1+\pi_t}, a_{i}) + \mathbb{E}_t \Theta_{t,t+1} V_{t+1}(\frac{p_{i}}{1+\pi_t}, a'_{i})
$$

$$
V^{A}_t(a_{i}) = \max_p \Pi(p, a_{i}) - \kappa + \mathbb{E}_t \Theta_{t,t+1} V_{t+1}(p, a'_{i})
$$

$$
V_t(p_{i},a_{i}) = \max \{V_t^{NA}(p_{i}, a_{i}), V_t^A(a_{i}) \}
$$

Where $pi_t = P_t/P_{t-1}$ is the inflation rate. If a firm doesnt adjust its prices its real price is deflated by the period's inflation rate.

Let $\chi(p_{i,t-1}, a_{it}) = \textbf{1}\{V^A(a_{it}) \geq V^{NA}(p_{i,t-1}, a_{it}) \}$ be the discrete decision on whether to change prices or not. 

Conditional on changing prices let the optimal pricing policy be given by $p'(a_{it})$. It only depends on $a_{it}$.


We can derive the aggregate Labour demand:

\begin{align}
L^d_t &= \int \frac{y_{it}}{e^{Z_t a_{it}}} di \\
&= \left(\int \frac{p_{it}^{-\epsilon}}{e^{Z_t a_{it}}} di \right) Y_t
\end{align}

Step 2 substitutes in the final demand for $y_{it}$ from the final goods aggregator.

## Cross Sectional Distribution

Let $\hat{\Omega}_t(p,a)$ be the cross sectional distribution at the beginning of $t$ after shocks have been realized but before any pricing decisions have been made. Let $\Omega_t(p,a)$ be the joint distribution of firms at the end of period $t$ i.e. after pricing decisions have been made.

# Consumers

Consumers are standard infinitely lived consumers who give us an euler equation and a static labour equation. We assume their utility functions are $u(C, L) =\log(C) - \zeta \frac{L_t^{1 + 1/\eta}}{1+1/\eta}$


The Euler equation is:

$$
\frac{1}{C_t} = \beta \frac{1 + r_t}{1+\pi_{t+1}} \frac{1}{C_{t+1}}
$$

The static labour optimality condition is

$$
L_t^{1/\eta} \zeta = \frac{w_t}{C_t}
$$

# Central bank

The central bank runs a taylor rule

$$
r_t = r_{t-1} + \phi_{\pi} (\pi_t - \pi^*) + \phi_{Y}(Y_t - Y^*)
$$

# Equilibrium

The recursive equilibrium is a set of functions $\{V^{A}, V^{NA}, p', \chi, C, Y,w, r, L, \pi, \Omega, \hat{\Omega}, Z\}$ satisfying the equations

\begin{align}
V^{NA}_{t-1}(p_i, a_i) &= \Pi(\frac{p_i}{1+\pi_{t-1}}, a_i) + \mathbb{E}_{t-1} \Theta_{t-1,t} V_{t}(\frac{p_i}{1+\pi_{t-1}}, a'_i) \\
V^{A}_{t-1}(a_i) &= \Pi(p'_{t-1}(a_i), a_i) - \kappa + \mathbb{E}_{t-1} \Theta_{t-1,t} V_{t}(p'_{t-1}(a_i), a'_i) \\
p'_{t-1}(a_i) &= \arg \max_p \Pi(p, a_i) + \mathbb{E}_{t-1} \Theta_{t-1,t} V_{t}(p, a'_i) \\
\chi_{t-1}(p_i, a_i) &= \textbf{1} \{V^A_{t-1}(a_i) > V^{NA}_{t-1}(p_i, a_i) \} \\
\hat{\Omega}_t(p_i, a'_i) &= \int \Omega_{t-1}(p_i, a_i) \Gamma(a'_i | a_i) d a_i \\
\Omega_t(p'_i, a_i) &= \int \chi(p_i, a_i) \textbf{1}\{p'(p_i, a_i) = p'_i \} \notag \\
&+ (1 - \chi(p_i, a_i)) \textbf{1}\{\frac{p_i}{1 + \pi_t} = p'_i\} \hat{\Omega}_t(p_i, a_i)dp_i \\
\frac{1}{C_t} &= \beta \frac{1 + r_t}{1+\pi_{t+1}} \frac{1}{C_{t+1}} \\
L_t^{1/\eta} \zeta &= \frac{w_t}{C_t} \\
r_t &= r_{t-1} + \phi_{\pi} (\pi_t - \pi^*) + \phi_{Y}(Y_t - Y^*) + \epsilon^r_{t} \\
1 &= \left(\int \int p_{i}^{1 - \epsilon} \Omega_t(p_i, a_i) dp_i da_i \right)^{\frac{1}{1-\epsilon}}  \\
Y_t &= C_t + \int \int \kappa \textbf{1}\{\chi_t(p_i, a_i) = 1\} \hat{\Omega}_t(p_i, a_i) dp_i da_i \\
L_t &= \int \int p_{i}^{-\epsilon} e^{-a_i} Y_t \Omega_t(p_i, a_i) dp_i da_i \\
Z_t &= \rho Z_{t-1} + \epsilon^z_{t}
\end{align}