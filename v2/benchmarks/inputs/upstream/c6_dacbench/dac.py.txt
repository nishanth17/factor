#!/usr/bin/env python3

def is_dac_starting_from_1(C):
  C = [int(m) for m in C]
  if len(C) == 0: return False
  if C[0] != 1: return False
  allowed = set([1])
  A = [1]
  for m in C[1:]:
    if any(m == 2*a or (m-a in allowed and 2*a-m in allowed) for a in A):
      allowed.add(m)
      A += [m]
      continue
    return False
  return True

def is_dasc_starting_from_1(C):
  C = [int(m) for m in C]
  if len(C) == 0: return False
  if C[0] != 1: return False
  allowed = set([1])
  A = [1]
  for m in C[1:]:
    if any(m == 2*a or (abs(m-a) in allowed and 2*a-m in allowed) for a in A):
      allowed.add(m)
      A += [m]
      continue
    return False
  return True

def test_misc():
  print('dac_misc')

  # miscellaneous short examples
  assert is_dac_starting_from_1([1])
  assert is_dac_starting_from_1([1,2])
  assert is_dac_starting_from_1([1,2,2])
  assert is_dac_starting_from_1([1,2,3])
  assert is_dac_starting_from_1([1,2,4])
  assert is_dac_starting_from_1([1,2,3,2])
  assert is_dac_starting_from_1([1,2,3,3])
  assert is_dac_starting_from_1([1,2,3,4])
  assert is_dac_starting_from_1([1,2,3,5])
  assert is_dac_starting_from_1([1,2,3,6])
  assert is_dac_starting_from_1([1,2,4,2])
  assert is_dac_starting_from_1([1,2,4,3])
  assert is_dac_starting_from_1([1,2,4,4])
  assert is_dac_starting_from_1([1,2,4,6])
  assert is_dac_starting_from_1([1,2,4,8])
  assert not is_dac_starting_from_1([])
  assert not is_dac_starting_from_1([0])
  assert not is_dac_starting_from_1([2])
  assert not is_dac_starting_from_1([1,1])
  assert not is_dac_starting_from_1([1,0])
  assert not is_dac_starting_from_1([1,3])
  assert not is_dac_starting_from_1([1,3,3])
  assert not is_dac_starting_from_1([1,2,1])
  assert not is_dac_starting_from_1([1,2,5])
  assert not is_dac_starting_from_1([1,2,3,1])
  assert not is_dac_starting_from_1([1,2,3,7])
  assert not is_dac_starting_from_1([1,2,3,8])
  assert not is_dac_starting_from_1([1,2,4,5])
  assert not is_dac_starting_from_1([1,2,4,7])

  # miscellaneous examples from 1992 montgomery
  assert is_dac_starting_from_1([1,2,3,4,5,9])
  assert is_dac_starting_from_1([1,2,3,4,6,7,13])
  assert is_dac_starting_from_1([1,2,3,4,7,8,15])
  assert is_dac_starting_from_1([1,2,3,4,5,8,9,17])
  assert is_dac_starting_from_1([1,2,3,6,9])
  assert is_dac_starting_from_1([1,2,3,5,8,13])
  assert is_dac_starting_from_1([1,2,3,5,10,15])
  assert is_dac_starting_from_1([1,2,3,6,9,15])
  assert is_dac_starting_from_1([1,2,3,4,7,10,17])
  assert is_dac_starting_from_1([1,2,3,5,6,11,17])
  assert is_dac_starting_from_1([1,2,3,5,7,10,17])
  assert is_dac_starting_from_1([1,2,3,5,7,12,17])
  assert not is_dac_starting_from_1([1,2,3,6,12,24,25,50,100,101])
  assert is_dac_starting_from_1([1,2,3,4,6,7,12,13,25,26,50,51,101])
  assert is_dac_starting_from_1([1,2,3,4,5,8,11,16,27,32,37,64,101]) # "right-to-left binary method"
  assert is_dac_starting_from_1([1,2,3,5,7,12,19,24,43,67,110,177,244,421])
  assert is_dac_starting_from_1([1,2,3,5,7,12,19,24,43,67,110,177,287,464,751])
  assert not is_dac_starting_from_1([1,2,3,4,7,10,11,9])
  assert is_dasc_starting_from_1([1,2,3,4,7,10,11,9])

def test_dac():
  import sys
  import random

  for r in range(1,20):
    print(f'dac {r}')
    sys.stdout.flush()
    for loop in range(100):
      C = [1]
      while len(C) < r:
        possible = set()
        for a in C:
          for b in C:
            if a == b or a-b in C:
              possible.add(a+b)
        C += [random.choice(list(possible))]
      assert is_dac_starting_from_1(C)
      assert is_dasc_starting_from_1(C)

def test_dasc():
  import sys
  import random

  for r in range(1,20):
    print(f'dasc {r}')
    sys.stdout.flush()
    for loop in range(100):
      C = [1]
      while len(C) < r:
        possible = set()
        for a in C:
          for b in C:
            if a == b or a-b in C:
              possible.add(a+b)
            if a+b in C:
              possible.add(abs(a-b))
        C += [random.choice(list(possible))]
      assert is_dasc_starting_from_1(C)

def test():
  test_misc()
  test_dac()
  test_dasc()

if __name__ == '__main__':
  test()
