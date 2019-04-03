# Question 1
We are provided the function
```
f(x) = sin^2(x) / x^2, x >= 0
```
In order to select x according to f we use importance sampling with
```
f(x) <= g(x) <= 1 (x < 1) + 1/x^2 (x >= 1)
```
The primitive of g is
```
G(x) = x (0 < x < 1) + (2 - 1/x) (x >= 1)
```
Note that the primitive spans (0,2). The inverse of G is
```
x = Ginv(u) = u (0 < u < 1) + 1/(2-u) (u >= 1), for u = 2R
```
We can now implement an RNG distributed according to f, see algorithm in q1.c
The algorithm also provides the integral over f for normalization.
We determine that Integral(f(x), x from 0 to infty) = 1.5707.
For comparison the analytic value is pi/2 = 1.5708.

The sampled points are shown in the histogram q1.pdf, where we have excluded points above 50.
The orange line is the true distribution for comparison.

#Q2
