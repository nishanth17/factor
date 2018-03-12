# factor

Fast prime factorization in Python. Factors most 50-60 digit numbers within a minute or so (with PyPy).  
The algorithm used depends on the size of the input

* `pollardPm1.py` contains an implementation of the large prime (two stage) variant of Pollard's _p-1_ algorithm.
* `pollardRho.py` contains an implementation of Pollard's _p-1_ algorithm with Brent's improvements. 
* `ecm.py` contains an implementation of Lenstra's elliptic curve factorization algorithm. It is inversionless (since it uses Montgomery coordinates), uses two stages, and uses Suyama's parametrization to generate elliptic curves. It also contains an implementation of Montgomery's PRAC algorithm for scalar multiplication (thanks Paul Zimmerman!) but this turned out to be slower than the usual double-and-add algorithm weirdly.
* `primeSieve.py` contains a bunch of prime sieves (Atkin, Eratosthenes, segmented Eratosthenes). Look at the [file](https://github.com/nishanth17/factor/blob/master/primeSieve.py) for specific benchmarks.

# Usage
All you have to do is run the file `factor.py`, enter a number, and hit Enter. Here's an example in terminal:

    python factor.py
    Enter a number: 15

    Factoring 15...
    Number of digits: 2
    Finding small prime factors...
    Prime factors found: 3, 5

    15 = 3^1 * 5^1

    Time: 5.00679016113e-05 s

# References
* A.O.L Atkin, D.J.Bernstein; [Prime Sieves using Binary Quadratic Forms](http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf); *Mathematics of Computation*, 73-246: 1023-30
* Peter L Montgomery; [Speeding the Pollard and Elliptical Methods of Factorization](http://modular.math.washington.edu/edu/124/misc/montgomery.pdf); *Mathematics of Computation* (Jan 1987), Issue 177: 243-264
* Montgomery, P.L.; [Evaluating Recurrences of the form X<sub>m+n</sub> = f(X<sub>m</sub>m, X<sub>n</sub>, X<sub>m-n</sub>) via Lucas Chains](http://cr.yp.to/bib/1992/montgomery-lucas.ps); Unpublished manuscript (Jan 1992)




    
