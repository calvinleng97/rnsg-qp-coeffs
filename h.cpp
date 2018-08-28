/** 
    h.cpp

    A Sequence of Quasipolynomials Arising From Random Numerical Semigroups

    Author: Calvin Leng

    Computes h_{n, d(n) - k} for sufficiently large n with respect to k, 
    i.e. the number of numerical semigroups S such that e(S) = d(n) - k,
    n not in S, and with minimal generating set A such that A < n / 2.
*/

#include <stdlib.h>
#include <iostream>
#include <vector>
#include <iomanip>
#include <set>
#include <ctime>
#include <cmath>
#include <climits>

using namespace std;

typedef unsigned long int ulli;
typedef set<int> SI;
typedef SI::iterator SetItr;
typedef vector<int> P;

/*
 * If `b` in binary is 1011 then `result` will be the
 * subset of `S` = {x_1, x_2, x_3, x_4} such that
 * `result` = {x_1, x_3, x_4}
 */
void generateSubset(const SI &S, const ulli b, SI &result)
{
  ulli c = 0;
  for (SetItr itr = S.begin(); itr != S.end(); itr++, c++)
    if ((b & (1 << c)) >> c)
      result.insert(*itr);
}

/*
 * Lower bound on I in Theorem 5.15
 *
 */
int p(const int k, const int b)
{
  return -2 * k - 1 + b;
}

int binom(const int n, int k)
{
  if (k > n) return 0;
  if (k * 2 > n) k = n - k;
  if (k == 0) return 1;

  int result = n;

  for (int i = 2; i <= k; ++i) {
    result *= (n - i + 1);
    result /= i;
  }

  return result;
}

bool in(int x, const SI &S)
{
  return S.find(x) != S.end();
}

/*
 * d(n) = size of X_n
 *
 */
int d(int n)
{
  return (n - 1) / 2 - n / 3;
}

/*
 * Checks if Theorem 5.5(iv) is satisfied
 * 
 */
bool isValidFixation(const vector<P> &pairs, const SI &fixation)
{
  for (unsigned int i = 0; i < pairs.size(); i++)
    if (!in(pairs[i][0], fixation) && !in(pairs[i][1], fixation))
      return false;

  return true;
}

void generateRemovingRange(const int b, const int m, SI &result)
{
  for (int i = 1; i <= b - 2 * m; i++)
    result.insert(i);
}

void setMinus(const SI &A, const SI &B, SI &result)
{
  for (SetItr itr = A.begin(); itr != A.end(); itr++) 
    if (!in(*itr, B)) 
      result.insert(*itr);
}

void updateCount(const SI &I, const int b, const int n, const int k, ulli &count)
{
  SI R;

  // Creates R = R \cup A(I) \cup B(I) \cup C(I)
  for (SetItr itr = I.begin(); itr != I.end(); itr++) {
    int a = *itr;

    // A(I)
    R.insert(b - 2*a);
    
    // B(I)
    if ( (b - a) % 2 == 0 ) R.insert((b - a) / 2); // B(I)

    // C(I)
    SetItr itr2 = itr;
    for (itr2++; itr2 != I.end(); itr2++) R.insert(b - a - *itr2);
  }

  int l = R.size() - I.size();

  if (l > k)
    return;

  // Creates P(I, R)
  vector<P> pairs;
  for (SetItr itr = I.begin(); itr != I.end(); itr++) {
    for (int x = 1; x < b - *itr - x; x++) {
      if (!(in(x, R) || in(b - *itr - x, R)))
      {
        P p;
        p.push_back(x);
        p.push_back(b - *itr - x);
        pairs.push_back(p);
      }
    }
  }

  // Creates R_c
  SI xn;
  generateRemovingRange(b, *(I.begin()), xn);
  SI R_c;
  setMinus(xn, R, R_c);

  // For all S in Powerset(R_c)
  for (ulli b1 = 0; b1 < pow(2, R_c.size()); b1++) {
    SI subset;
    generateSubset(R_c, b1, subset);
    int l_ = l + subset.size();
     
    // Condition for Theorem 5.5(iv) 
    if (l_ <= k && isValidFixation(pairs, subset))
      count += binom(d(n) - xn.size(), k - l_);
  }
}

int main(int argc, char **argv)
{
  if (argc < 3) {
    cout << "Error: must supply at least 2 arguments: n and k." << endl;
    return 1;
  }

  int n = atoi(argv[1]);
  int k = atoi(argv[2]);
  int b = n % 3;

  if (n <= 24 * k + 12 - 8 * b) {
    cout << "Error: n must strictly greater than " 
         << 24 * k + 12 - 8 * b << " for k = " << k << endl;
    return 1;
  }

  clock_t start_t, end_t;
  start_t = clock();
  SI insertingRange;
  ulli count = 0;

  for (int i = 0 - (n % 3 == 0); i >= p(k, b); i--)
    insertingRange.insert(i);

  // For I in Powerset({p_n(k), p_n(k) + 1, ..., 0})
  for (ulli b1 = 0; b1 < pow(2, insertingRange.size()); b1++)
  {
    SI I;
    generateSubset(insertingRange, b1, I);
    updateCount(I, b, n, k, count);
  }

  end_t = clock();
  cout << "There are " << count << " numerical semigroups of embedding dimension " << d(n) - k 
       << " with minimal generating set bounded above by " << n / 2 << " such that " << n << " is not in the semigroup, i.e.\n" 
       << "h_{" << n << ", " << d(n) - k << "} = " << count << endl;
  cout << (float) ((float) end_t - (float) start_t) << "ms to run." << endl;

  return 0;
}
