# Parton Showers
We consider the Parton shower with the processes
`q->qg` and `g->gg`.
These have the rate pabctions
```
fabc(p,z) = alpha_s / (2*pi) Pabc(z) / p^2
```
where
`Pqqg(z) = 4/3*(1+z^2)/(1-z)` and `Pggg(z)=3(1-z(1-z))^2/(z*(1-z))`.

We note that both `Pqqg` and `Pggg` are maximized for the limit `z=1-1/E_a`
and use veto sampling with
```
gabc(p,z) = alpha_s / (2*pi) Pabc(1-1/E_a) / p^2
```
The veto algorithm is implemented in the file `parton.c` as

```C
void veto(double ea, double pmax, double (*pabc)(double z), gsl_rng* r, double *pout, double *zout)
{
        double term, z, y, ft, p, p2;

        // determine time
        y = 1;
        ft = 0;
        p2 = pmax*pmax;
        do {
                term = GINT(p2,pabc)-log(gsl_rng_uniform_pos(r));
                p2 = GINV(term,pabc);
                p = sqrt(p2);
                if (2*p < ea)
                {
                        z = p/ea + gsl_rng_uniform_pos(r)*(1.-2.*p/ea);
                        y = gsl_rng_uniform_pos(r)*GF(p2,pabc);
                        ft = 0.5*M_1_PI*ALPHA/p2*pabc(z);
                }
        } while (y > ft);
        *pout = p;
        *zout = z;

        return;
}
```

The shower is handled by calling the veto algorithm recursively and counting
the number of `q->qg` events as well as the number of added gluons.
Sampling 1000000 we obtain
|`E_0`		| `avg(q->qg)`	| `sd(q->qg)`	| `avg(#g)`	| `sd(#g)`|
|-----------|---------------|---------------|-----------|---------|
|;50|2.03764|0.928995|4.8761|4.1521|
|100|2.40562|1.06838|7.58461|6.28167|
|200|2.8382|1.21279|12.0711|9.57787|
|400|3.33221|1.35779|19.5401|14.8033|

