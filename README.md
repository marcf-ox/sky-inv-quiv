<section>
  <h2>Summary</h2>
   This package implements in Python an algorithm to compute the skyscaper invariant [2] over a finite grid poset. This algorithm relies on the ideas described in [1] and works for coefficients in any field. In particular, when the coefficients are in
- $\mathbb Q$ — we use the cfraction library for exact computations
- $\mathbb R$ or $\mathbb C$ — we use NumPy/SciPy, which may create errors for matrices that are not well conditioned
- $\mathbb Z/2\mathbb Z$ — we recommend using instead the algorithm at https://github.com/JanJend/Skyscraper-Invariant
</section>

<section>
  <h2>
    Installation
  </h2>
    It may be installed via

```
pip install git+https://github.com/marcf-ox/sky-inv-quiv
```

This requires Python >=3.9.0, Numpy>=2.1.1, Scipy>=1.9.0, cfractions>=2.3.1, Networkx>=2.5.1, Matplotlib>=3.9.1.
</section>

<section>
  <h2>
    Example
  </h2> 

 A Jupyter tutorial notebook is available at:

<a href = https://github.com/marcf-ox/HNcode/notebook/tuto.ipynb> Tutorial notebook </a>


</section> 


<section>
  <h2>References</h2>
  <ol>
    <li>
      Cheng, Chi-Yu. (2024). <cite>A deterministic algorithm for Harder–Narasimhan filtrations for representations of acyclic quivers</cite>. <em>Algebra & Number Theory</em>, <strong>18</strong>(2), 319–347. 
      <a href="https://doi.org/10.2140/ant.2024.18.319" target="_blank" rel="noopener">https://doi.org/10.2140/ant.2024.18.319</a>
    </li>
    <li>
      Fersztand, Marc, Jacquard, Emile, Nanda, Vidit, & Tillmann, Ulrike. (2024). <cite>Harder–Narasimhan filtrations of persistence modules</cite>. <em>Transactions of the London Mathematical Society</em>, <strong>11</strong>(1), e70003. 
      <a href="https://doi.org/10.1112/tlm3.70003" target="_blank" rel="noopener">https://doi.org/10.1112/tlm3.70003</a>
    </li>
  </ol>
</section>
