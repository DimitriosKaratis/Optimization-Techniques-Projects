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

**Results & Observations:**  
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

**Results & Observations:**  
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

**Results & Observations:**  
- For convergence, the step size **γk** must remain within strict theoretical bounds depending on the Lipschitz constant of the gradient.  
- Too small a step size resulted in extremely slow convergence, while too large a step size caused divergence.  
- Experimental results matched the theoretical predictions closely: oscillations and divergence appeared exactly when γk exceeded the upper bound.  
- The analysis highlighted the **critical role of adaptive step-size selection** in gradient-based optimization.  

---

### 🧬 **Project 4 – Genetic Algorithms (Final Project)**  
The final project focused on **evolutionary optimization techniques**, specifically the design and implementation of a **Genetic Algorithm (GA)** to solve challenging optimization problems.  

#### 📌 Tasks  
- **Problem Definition**: Apply GA to a non-trivial optimization problem (function minimization or classification), where classical deterministic methods may fail due to non-convexity or multiple local minima.  
- **Chromosome Encoding**: Decide how solutions are represented (binary, real-valued, or custom structures depending on the problem).  
- **Population Initialization**: Generate diverse initial candidate solutions to ensure sufficient exploration of the search space.  
- **Genetic Operators**:  
  - **Selection** – Implement parent selection mechanisms (e.g., roulette wheel, tournament selection).  
  - **Crossover** – Combine parents to produce new offspring (single-point, two-point, or uniform crossover).  
  - **Mutation** – Introduce random modifications to preserve genetic diversity and avoid premature convergence.  
- **Parameter Tuning**: Study how parameters such as population size, crossover probability, and mutation rate influence the GA’s convergence and solution quality.  
- **Stopping Criteria**: Define when the algorithm should terminate (fixed iterations, fitness threshold, or stagnation).  
- **Comparative Analysis**: Compare the performance of GA with classical methods studied in earlier projects (Steepest Descent, Newton, Levenberg–Marquardt).  

#### 📊 Results & Observations  
- The GA successfully solved optimization problems where deterministic methods often became trapped in local minima.  
- A **balanced mutation rate** was essential:  
  - Too low → lack of diversity, risk of premature convergence.  
  - Too high → random search behavior, slow convergence.  
- The **crossover probability** strongly affected exploitation of good candidate solutions.  
- Larger populations generally improved robustness but increased computational cost.  
- Compared to gradient-based methods, the GA required more iterations but achieved **better global optima** in multi-modal landscapes.  
- The results highlighted the **strength of evolutionary computation** for global optimization, especially when problem structure is complex or derivative information is unavailable.  

---

## ⚙️ Tools & Implementation
- **MATLAB** was used for all implementations and simulations.  
- Detailed comments and explanations are provided in each `.m` file.  
- Reports are written in Greek, while the README provides an English overview.  
