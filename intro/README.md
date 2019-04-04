# Question 1
We are provided the function
```
f(x) = sin^2(x) / x^2, x >= 0
```
In order to select x according to `f` we use importance sampling with
```
f(x) <= g(x) <= 1 (x < 1) + 1/x^2 (x >= 1)
```
The primitive of `g` is
```
G(x) = x (0 < x < 1) + (2 - 1/x) (x >= 1)
```
Note that the primitive spans `(0,2)`. The inverse of `G` is
```
x = Ginv(u) = u (0 < u < 1) + 1/(2-u) (u >= 1), for u = 2R
```
We can now implement an RNG distributed according to `f`, see algorithm `rng_sinc` in `q1.c`.
The algorithm also provides the integral over `f` for normalization.
We determine that the numerical integral is
```
NInt(f(x), x: 0 --> infty) = 1.5707
```
For comparison the analytic value is `pi/2 = 1.5708`.

The sampled points are shown in the histogram `q1.pdf`,
where we have excluded points above 50.
The orange line is the true distribution for comparison.

# QUESTION 2
The next function is
```
f(x) = cos(pi*x/2)
```
## 2.a
The analytic integral of `f` is
```
Int(f(x), x: 0 --> 1) = 2/pi = 0.637
```
## 2.b
We use hit-and-miss with flat sampling
```
f(x) <= 1
```
See function `mci_cosflat` in `q2.c` for implementation. The convergence along the
sampling of a million points is visible in `q2.pdf` (blue line).
From the log-log plot we see that the convergence is roughly a power-law.
We find
```
NInt(f(x), x: 0 --> 1) = 0.6368
```

## 2.c
We now use importance sampling with
```
f(x) <= g(x) = 1- x^2, 0 <= x <= 1
```
The primitive of `g` is
```
G(x) = x - 1/3 x^3
```
with an inverse provided in the exercise. Note that `G(x)` spans `(0,2/3)`.
From importance sampling the computed integral value is
```
NInt(f(x), x: 0 --> 1) = 0.6364
```
The convergence is visible in `q2.pdf` (orange line).
Due to the fact that importance sampling, by design, should sample more
in the important areas of the function we would expect this to converge faster.
However, in this specific case, a faster trend is not visible in the obtained
estimates.
