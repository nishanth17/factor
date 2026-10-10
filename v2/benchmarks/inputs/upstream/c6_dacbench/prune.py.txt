#!/usr/bin/env python3

import ladder
import tally

def search(target, target_len, too_small):
    ch = [1,2,3]
    chain = []
    chainlen = 0 
    while ch[-1] != target:
        if chainlen > target_len or ch[-1] > target or ch[-1] <= too_small[chainlen+2]: #over, over or undershoot, thus must be pruned
            if chainlen == 0:
                return None, target_len
            last_bit = chain.pop()
            ch.pop()
            remind = ch[-1] - ch[-2]
            chainlen -= 1
            count = 0
            while last_bit == 1: #tried both children, keep moving up until we find a 0, then try right
                    count += 1
                    if chainlen == 0:
                        return None, target_len
                    chainlen -= 1
                    last_bit = chain.pop()
                    ch.pop()
                    remind = ch[-1] - ch[-2]
            if last_bit == 0: #tried left child + tree, and try right 
                tally.node('right')
                ch.append(remind + ch[-1])
                chain.append(1)
                chainlen += 1
        else: #Take left hand node
            tally.node('left')
            ch.append(ch[-1] + ch[-2])
            chain.append(0)
            chainlen += 1
    else: 
        return chain, chainlen
    
def chain(n):
    if n < 3: return ladder.chain(n)
    best = None
    target_len = -1
    while best is None:
        target_len += 1
        f0 = 1
        f1 = 1
        too_small = [0]*(target_len + 3)
        for j in range(target_len + 2, -1, -1):
            too_small[j] = int((n-1)/f1)
            f2 = f0 + f1
            f0 = f1
            f1 = f2
        best,bestlen = search(n,target_len,too_small)
    C = [1,2,3]
    r0,r1,r2 = C
    for bestj in best:
      if bestj: r0,r1 = r1,r0
      r0,r1,r2 = r1,r2,r1+r2
      C += [r2]
    return C

def test_primes():
  import sys
  from primes import primes

  for lim in 10,100,1000,10000:
    P = primes(lim)
    tally.clear()
    total = sum(len(chain(n))-1 for n in P)
    print(f'prune primes limit {lim} primes {len(P)} primesum {sum(P)} chainlen {total} searchnodes {tally.nodes}')
    sys.stdout.flush()

def test_dac():
  import sys
  import dac

  assert chain(29) == [1,2,3,5,8,13,21,29]

  for n in range(1,1025):
    C = chain(n)
    # print(f'prune dac {n} {C}')
    # sys.stdout.flush()
    assert C[-1] == n
    assert dac.is_dac_starting_from_1(C)
    if n&(n-1) == 0:
      print(f'prune dac {n}')
      sys.stdout.flush()

def test():
  test_dac()
  test_primes()

if __name__ == '__main__':
  test()
