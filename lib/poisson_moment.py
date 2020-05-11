'''
Contains:
	- touchard: function to compute moment of Poisson distribution 
'''
def stirling2 ( n, m ):

#*****************************************************************************80
#
## STIRLING2 computes the Stirling numbers of the second kind.
#
#  Discussion:
#
#    S2(N,M) represents the number of distinct partitions of N elements
#    into M nonempty sets.  For a fixed N, the sum of the Stirling
#    numbers S2(N,M) is represented by B(N), called "Bell's number",
#    and represents the number of distinct partitions of N elements.
#
#    For example, with 4 objects, there are:
#
#    1 partition into 1 set:
#
#      (A,B,C,D)
#
#    7 partitions into 2 sets:
#
#      (A,B,C) (D)
#      (A,B,D) (C)
#      (A,C,D) (B)
#      (A) (B,C,D)
#      (A,B) (C,D)
#      (A,C) (B,D)
#      (A,D) (B,C)
#
#    6 partitions into 3 sets:
#
#      (A,B) (C) (D)
#      (A) (B,C) (D)
#      (A) (B) (C,D)
#      (A,C) (B) (D)
#      (A,D) (B) (C)
#      (A) (B,D) (C)
#
#    1 partition into 4 sets:
#
#      (A) (B) (C) (D)
#
#    So S2(4,1) = 1, S2(4,2) = 7, S2(4,3) = 6, S2(4,4) = 1, and B(4) = 15.
#
#
#  First terms:
#
#    N/M: 1    2    3    4    5    6    7    8
#
#    1    1    0    0    0    0    0    0    0
#    2    1    1    0    0    0    0    0    0
#    3    1    3    1    0    0    0    0    0
#    4    1    7    6    1    0    0    0    0
#    5    1   15   25   10    1    0    0    0
#    6    1   31   90   65   15    1    0    0
#    7    1   63  301  350  140   21    1    0
#    8    1  127  966 1701 1050  266   28    1
#
#  Recursion:
#
#    S2(N,1) = 1 for all N.
#    S2(I,I) = 1 for all I.
#    S2(I,J) = 0 if I < J.
#
#    S2(N,M) = M * S2(N-1,M) + S2(N-1,M-1)
#
#  Properties:
#
#    sum ( 1 <= K <= M ) S2(I,K) * S1(K,J) = Delta(I,J)
#
#    X^N = sum ( 0 <= K <= N ) S2(N,K) X_K
#    where X_K is the falling factorial function.
#
#  Licensing:
#
#    This code is distributed under the GNU LGPL license.
#
#  Modified:
#
#    26 February 2015
#
#  Author:
#
#    John Burkardt
#
#  Parameters:
#
#    Input, integer N, the number of rows of the table.
#
#    Input, integer M, the number of columns of the table.
#
#    Output, integer S2(N,M), the Stirling numbers of the second kind.
#
  import numpy as np

  s2 = np.zeros ( ( n, m ) ).astype(int) 

  if ( n <= 0 ):
    return s2

  if ( m <= 0 ):
    return s2

  s2[0,0] = 1
  for j in range ( 1, m ):
    s2[0,j] = 0

  for i in range ( 1, n ):

    s2[i,0] = 1

    for j in range ( 1, m ):
      s2[i,j] = int(( j + 1 ) * s2[i-1,j] + s2[i-1,j-1])

  return s2

def touchard(n,x):
    '''
    touchard(n,x) gives the nth-moment of the Poisson distribution centered in x
    '''
    if n < 0:
        return NULL
    if n == 0:
        return 1
    stirling_tb=stirling2(n, n)
    return sum(stirling_tb[n-1,k] * x ** (k+1) for k in range(n))

