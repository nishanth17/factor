#!/usr/bin/env python3

import tally

def chain(n):
  n = int(n)
  assert n >= 1
  result = set([n])
  two = 2
  while two <= n:
    tally.node('step')
    result.add(n // two)
    result.add((n+two-1) // two)
    two *= 2
  return sorted(result)

def test_misc():
  print('ladder misc')

  assert chain(1) == [1]
  assert chain(2) == [1,2]
  assert chain(3) == [1,2,3]
  assert chain(4) == [1,2,4]
  assert chain(5) == [1,2,3,5]
  assert chain(6) == [1,2,3,6]
  assert chain(29) == [1,2,3,4,7,8,14,15,29]
  assert chain(101) == [1,2,3,4,6,7,12,13,25,26,50,51,101] # example from 1992 montgomery

def test_lengths():
  import dac

  print('ladder lengths')
  
  def floorlog2(n):
    result = 0
    while (2<<result) <= n:
      result += 1
    return result

  def predictedlen(n):
    # formulas from 1992 montgomery, section 2
    if n == 1: return 0
    if n%2 == 0: return 1+predictedlen(n//2)
    return floorlog2(n)+floorlog2((2*n)//3)

  for n in range(1,10000):
    C = chain(n)
    assert C[-1] == n
    assert dac.is_dac_starting_from_1(C)
    assert len(C)-1 == predictedlen(n)

def test_primes():
  from primes import primes
  import sys

  for lim in 10,100,1000,10000,100000,1000000:
    P = primes(lim)
    tally.clear()
    total = sum(len(chain(n))-1 for n in P)
    print(f'ladder primes limit {lim} primes {len(P)} primesum {sum(P)} chainlen {total} searchnodes {tally.nodes}')
    sys.stdout.flush()

def test():
  test_misc()
  test_lengths()
  test_primes()

if __name__ == '__main__':
  test()
