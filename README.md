# Optimization Techniques Projects

This repository contains the implementation of multiple projects in **Optimization Techniques**, focusing on classical deterministic methods, gradient-based approaches, convergence analysis, and evolutionary optimization (Genetic Algorithms).  

The work was developed as part of the *Optimization Techniques* course at the **Department of Electrical and Computer Engineering, Aristotle University of Thessaloniki**.  

---

### 🔹 **Project 1 – Unimodal Function Optimization**
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
Application of **Genetic Algorithms (GA)** to complex optimization problems:  
- Chromosome representation and initialization strategies.  
- Selection, crossover, and mutation operators.  
- Convergence analysis and parameter tuning (population size, mutation rate, crossover probability).  
- Comparison of GA performance with classical methods from previous projects.  

**Results & Observations:**  
- GAs were effective in solving non-convex and multimodal problems where deterministic methods failed or got trapped in local minima.  
- Proper tuning of **population size** and **mutation rate** was crucial: too low mutation reduced exploration, while too high mutation prevented convergence.  
- Crossover probability strongly influenced the algorithm’s ability to exploit good solutions.  
- Compared to classical methods (Steepest Descent, Newton, etc.), GAs required more iterations but achieved **better global solutions** in complex landscapes.  
- The results confirmed the strength of **evolutionary methods in global optimization**, despite their slower convergence compared to deterministic approaches.  

---

## ⚙️ Tools & Implementation
- **MATLAB** was used for all implementations and simulations.  
- Detailed comments and explanations are provided in each `.m` file.  
- Reports are written in Greek, while the README provides an English overview.  
