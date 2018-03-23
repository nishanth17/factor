# factor

Fast prime factorization in Python. Factors most 50-60 digit numbers within a minute or so (with PyPy).  
The algorithm used depends on the size of the input

* `pollardPm1.py` contains an implementation of the large prime (two stage) variant of Pollard's _p-1_ algorithm.
* `pollardRho.py` contains an implementation of Pollard's Rho algorithm with Brent's improvements. 
* `ecm.py` contains an implementation of Lenstra's elliptic curve factorization algorithm. It is inversionless (since it uses Montgomery coordinates), uses two stages, and uses Suyama's parametrization to generate random elliptic curves. It also contains an implementation of Montgomery's PRAC algorithm for scalar multiplication (thanks Paul Zimmerman!) but this turned out to be slower than the usual double-and-add algorithm weirdly.
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

and another...

	Enter number: 37897387397398739739826929827929827927927762729872987928

	Factoring 37897387397398739739826929827929827927927762729872987928...
	Number of digits: 56
	Finding small prime factors...
	Prime factors found: 2, 3
	Factoring 1579057808224947489159455409497076163663656780411374497 with ECM...
	Number of digits: 55
	Bounds: 250000 128992510
	Sieving primes...
	Stage 2 found factor!
	Found factor 67246307
	Factoring 67246307...
	Number of digits: 8
	67246307 is prime!
	Factoring 23481702991138940747474138758238071923617408171...
	Number of digits: 47
	Factoring 23481702991138940747474138758238071923617408171 with ECM...
	Number of digits: 47
	Bounds: 50000 12746592
	Sieving primes...
	Tried 40 random curves...
	Tried 80 random curves...
	Tried 120 random curves...
	Tried 160 random curves...
	Stage 2 found factor!
	Found factor 4788272261623351
	Factoring 4788272261623351...
	Number of digits: 16
	4788272261623351 is prime!
	Factoring 4904003303934522319753958187821...
	Number of digits: 31
	4904003303934522319753958187821 is prime!

	37897387397398739739826929827929827927927762729872987928 = 2^3 * 3^1 * 67246307^1 * 4788272261623351^1 * 4904003303934522319753958187821^1

	Time: 24.7774269581 s

# References
* A.O.L Atkin, D.J.Bernstein; [Prime Sieves using Binary Quadratic Forms](http://www.ams.org/journals/mcom/2004-73-246/S0025-5718-03-01501-1/S0025-5718-03-01501-1.pdf); *Mathematics of Computation*, 73-246: 1023-30
* Peter L Montgomery; [Speeding the Pollard and Elliptical Methods of Factorization](http://modular.math.washington.edu/edu/124/misc/montgomery.pdf); *Mathematics of Computation* (Jan 1987), Issue 177: 243-264
* Montgomery, P.L.; [Evaluating Recurrences of the form <i>X<sub>m+n</sub></i> = <i>f(X<sub>m</sub>, X<sub>n</sub>, X<sub>m-n</sub>)</i> via Lucas Chains](http://cr.yp.to/bib/1992/montgomery-lucas.ps); Unpublished manuscript (Jan 1992)




    
