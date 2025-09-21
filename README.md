# Optimization Techniques Projects

This repository contains the implementation of multiple projects in **Optimization Techniques**, focusing on classical deterministic methods, gradient-based approaches, convergence analysis, and evolutionary optimization (Genetic Algorithms).  

The work was developed as part of the *Optimization Techniques* course at the **Department of Electrical and Computer Engineering, Aristotle University of Thessaloniki**.  

---

### 🔹 **Project 1 – Optimization of Single-Variable Functions**
Implementation and comparison of classical 1D optimization algorithms:
- **Dichotomy Method**  
- **Golden Section Method**  
- **Fibonacci Method**  
- **Dichotomy with Derivatives**  

Each method was applied to three different test functions, analyzing:
- Number of function evaluations vs. accuracy parameters.  
- Interval reduction across iterations.  
- Comparative efficiency of the methods.  

#### 📊 Results & Observations  
- The **Golden Section** and **Fibonacci** methods were the most efficient in terms of convergence speed and robustness.  
- The **Dichotomy method** was simpler but required more function evaluations due to repeated interval splitting.  
- When derivative information was available, the **Dichotomy with Derivatives** method accelerated convergence, but its performance depended heavily on numerical stability of derivative estimation.  
- Overall, Fibonacci offered the best balance between accuracy and computational cost.  

---

### 🔹 **Project 2 – Multivariable Optimization Methods**
Study of optimization in 2D functions using different algorithms and step-size strategies:  
- **Steepest Descent Method**  
- **Newton’s Method**  
- **Levenberg–Marquardt Method**  

Step-size selection strategies explored:
- Fixed step size.  
- Line search with **Golden Section**.  
- **Armijo Rule**.  

#### 📊 Results & Observations  
- **Steepest Descent** converged slowly and was highly sensitive to the choice of step size; Golden Section and Armijo improved stability compared to a fixed step.  
- **Newton’s Method** provided fast convergence when the Hessian was well-conditioned, but failed or diverged when the Hessian was ill-conditioned.  
- **Levenberg–Marquardt** successfully combined the speed of Newton with the stability of gradient descent, especially in difficult regions.  
- Among the step-size strategies, **Armijo rule** provided the best trade-off between convergence speed and robustness.  

---

### 🔹 **Project 3 – Convergence of Gradient Methods**
In-depth analysis of the **convergence properties of the Steepest Descent method**, focusing on:  
- The influence of step size (γk).  
- Conditions for convergence and divergence.  
- Mathematical justification for the step-size bounds.  
- Empirical validation through simulations with different initial conditions.  

#### 📊 Results & Observations   
- For convergence, the step size **γk** must remain within strict theoretical bounds depending on the Lipschitz constant of the gradient.  
- Too small a step size resulted in extremely slow convergence, while too large a step size caused divergence.  
- Experimental results matched the theoretical predictions closely: oscillations and divergence appeared exactly when γk exceeded the upper bound.  
- The analysis highlighted the **critical role of adaptive step-size selection** in gradient-based optimization.  

---

### 🧬 **Project 4 – Genetic Algorithms (Traffic Network Optimization)**  

The final project applied **evolutionary optimization techniques** to a real-world inspired problem:  
**minimizing total travel time in a congested traffic network.**  

#### 🚦 Problem Description  
- The network is modeled as a directed graph where **nodes** are intersections and **edges** are roads.  
- Each road *i* has:  
  - **Baseline travel time**: *t₁, t₂, …* (minutes, under negligible traffic).  
  - **Capacity**: *c₁, c₂, …* (vehicles per minute, maximum flow).  
  - **Flow variable**: *x₁, x₂, …* (vehicles per minute on road *i*).  
  - **Congestion factor**: *a₁, a₂, …*, a constant depending on road type.  

$$
T_i(x_i) = t_i + \frac{a_i x_i}{1 - \frac{x_i}{c_i}} \quad \text{[minutes]}
$$

- Properties:  
  - As $x_i \to 0$, $T_i(x_i) \to t_i$ (light traffic).  
  - As $x_i \to c_i$, $T_i(x_i) \to +\infty$ (congestion → gridlock). 

- **Constraints**:  
  - The total incoming flow at each intersection equals the total outgoing flow (no accumulation or deficits).  
  - The total input flow for the network is fixed at \(V = 100\) vehicles/minute.  

- **Objective**:  
  $$
  \text{Minimize} \quad T_\text{total} = \sum_{i=1}^{n} T_i(x_i)
  $$ 
  subject to flow conservation and capacity constraints.  

#### 📌 Tasks  
- **Mathematical formulation** of the traffic optimization problem.  
- **Chromosome representation** of road flows $(x_1, x_2, \dots, x_n)$.  
- **Population initialization** ensuring feasibility (flows respect conservation laws and capacities).  
- **Genetic Operators**:  
  - **Selection** – roulette wheel and tournament selection.  
  - **Crossover** – combining road-flow distributions between parent solutions.  
  - **Mutation** – introducing small random changes to flows for exploration.  
- **Parameter tuning**: analyze the impact of population size, mutation rate, and crossover probability.  
- **Comparative study**: evaluate GA performance against classical methods from earlier projects (Steepest Descent, Newton, Levenberg–Marquardt).  

#### 📊 Results & Observations  
- The **GA successfully minimized total travel time** across the network while respecting road capacities and flow conservation.  
- **Mutation rate** was crucial:  
  - Too low → premature convergence to suboptimal solutions.  
  - Too high → excessive randomness, poor convergence.  
- **Crossover probability** balanced exploration vs. exploitation; optimal values helped combine partial solutions efficiently.  
- Larger populations improved solution quality but increased runtime.  
- Compared to deterministic methods, GA required more iterations but avoided local minima and achieved **better global traffic flow distribution**.  
- The project demonstrated the **practical advantage of evolutionary computation** for large-scale, non-convex optimization problems like traffic networks.  

---

## ⚙️ Tools & Implementation
- **MATLAB** was used for all implementations and simulations.  
- Detailed comments and explanations are provided in each `.m` file.  
- Reports are written in Greek, while the README provides an English overview.  
