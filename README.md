# Weyl growth and tracial limits for arithmetic translations

R. Arun Chandru. Manuscript dated 27 September 2026.

## Abstract

We determine the growing-window tracial pressure of the actual prime-power
translation operator attached to a fixed holomorphic cuspidal newform of
trivial character. Every intermediate translation is killed at the interval
boundary. At fixed coupling s the logarithmic trace is
s^2 ell^2/12 + O(ell sqrt(log ell)), with a bounded lower error; the
corresponding fixed-window Weyl coefficient has quadratic logarithmic
growth ell^2/48. A uniform coupling estimate gives a large-deviation law at
every diverging intermediate scale up to ell^2. At the central scale ell,
the arithmetic trace converges to the square-root variance-profile Toeplitz
law, with moments 1/3 and 4/15 and logarithmic tail coefficient -3 in both
directions. An exact weighted autocorrelation identity explains the pressure
constant. We prove global L2 stability modulo affine phase, with optimal
square-root exponent and sharp small-deficit squared-distance coefficient 6.
An explicit taper family preserves total leading prime variance while
changing the boundary pressure. Higher prime powers and ramification are
restored at bounded logarithmic cost, and fixed tuples of inequivalent forms
give an orthogonally invariant, noncommuting colored limit.

## Contents

- manuscript.pdf: complete manuscript, including references.
- manuscript.tex: single, self-contained LaTeX source; bibliography embedded.
- verify_identities.py: optional exact finite algebra checks.
- README.md: this file.


## Scope and attribution

The fixed-window Weyl formula, Rankin--Selberg prime number theorem,
Deligne bound, general Toeplitz limit theory, and Gaussian concentration
are cited inputs. The manuscript distinguishes them from the actual
arithmetic growing-window comparison, weighted profile geometry and
stability, profile-specific tail analysis, and taper results developed here.

Forms, tuples, and coupling bounds are fixed in the stated uniform
estimates. There is no assertion of uniformity in conductor or weight, nor
of a window-uniform high-energy Weyl onset. The trace is the specified
normalized translation-algebra trace, not ordinary Hilbert-space heat.

The cosine parameter interval is exact for the sufficient Fourier-sign
certificate used in the proof, not asserted to be the maximal interval
of flat-profile optimality. The global stability constant is finite and
absolute but not given numerically; 6 is the sharp small-deficit
squared-distance coefficient, not a claimed sharp global constant.

## Compilation

Use a current LaTeX distribution providing amsart, AMS mathematics,
mathtools, mathrsfs, lmodern, geometry, microtype, xurl, and hyperref.
No shell escape or BibTeX run is needed.

From this folder:

~~~sh
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
pdflatex -interaction=nonstopmode -halt-on-error manuscript.tex
~~~

If cross-references still request another run, repeat the command.
Alternatively:

~~~sh
tectonic manuscript.tex
~~~


## Optional exact checks

Python 3 with the standard library is sufficient:

~~~sh
python verify_identities.py
~~~

The script checks rational moment integrals and normalizations, discrete
overlap identities, polynomial instances of the transverse quadratic
formulas, and the algebra of the taper certificate. 

These finite checks do not establish the arithmetic prime number theorem,
infinite-dimensional stability, large-deviation limits, or other
asymptotic assertions. Their proofs and cited inputs are in the manuscript;
no limiting theorem rests on numerical evidence.
