# R estimation code for "Incorporating Time-Varying Covariates in a Simple Mixture Model for Discrete-Time Duration Data"

Pete Fader and Bruce Hardie provide a technical note https://brucehardie.com/notes/037/ on incorporating time-varying covariates for a discrete-time duration mixture model. It is an extention to the (shifted) Beta-Geometric Distribution, now with covariates. They call it **G2G+covariates model**.

However, the code on how to estimate this model is not available.  I created the MLE estimation code using R.

## How to use

The user needs to prepare the data in _person-period_ format.  For example:

```
> print(Scania_PersonPeriod_Train[1:50,], n=50)
# A tibble: 50 × 13
      id enter  exit event birthdate sex    parish ses   immigrant spell   age  year foodprices
   <int> <dbl> <int> <dbl>     <dbl> <fct>  <fct>  <fct> <fct>     <dbl> <dbl> <dbl>      <dbl>
 1   415     0     1     0      1805 female 3      lower yes          28    50  1855      0.284
 2   415     1     2     0      1805 female 3      lower yes          28    51  1856      0.091
 3   415     2     3     0      1805 female 3      lower yes          28    52  1857      0.084
 4   415     3     4     0      1805 female 3      lower yes          28    53  1858     -0.145
 5   415     4     5     0      1805 female 3      lower yes          28    54  1859     -0.241
 6   415     5     6     0      1805 female 3      lower yes          28    55  1860     -0.013
 7   415     6     7     0      1805 female 3      lower yes          28    56  1861      0.102
 8   415     7     8     0      1805 female 3      lower yes          28    57  1862     -0.052
 9   415     8     9     0      1805 female 3      lower yes          28    58  1863     -0.122
10   415     9    10     0      1805 female 3      lower yes          28    59  1864     -0.103
11   415    10    11     0      1805 female 3      lower yes          28    60  1865      0.12 
12   415    11    12     0      1805 female 3      lower yes          28    61  1866     -0.086
13   415    12    13     0      1805 female 3      lower yes          28    62  1867      0.099
14   415    13    14     0      1805 female 3      lower yes          28    63  1868      0.299
15   415    14    15     0      1805 female 3      lower yes          28    64  1869     -0.171
16   415    15    16     0      1805 female 3      lower yes          28    65  1870     -0.151
17   415    16    17     0      1805 female 3      lower yes          28    66  1871     -0.093
18   415    17    18     0      1805 female 3      lower yes          28    67  1872     -0.074
19   415    18    19     0      1805 female 3      lower yes          28    68  1873      0.204
20   415    19    20     0      1805 female 3      lower yes          28    69  1874      0.194
21   415    20    21     0      1805 female 3      lower yes          28    70  1875      0.003
22   415    21    22     0      1805 female 3      lower yes          28    71  1876      0.021
23   415    22    23     0      1805 female 3      lower yes          28    72  1877      0.079
24   415    23    24     0      1805 female 3      lower yes          28    73  1878     -0.152
25   415    24    25     0      1805 female 3      lower yes          28    74  1879     -0.167
26   415    25    26     0      1805 female 3      lower yes          28    75  1880      0.221
27   415    26    27     0      1805 female 3      lower yes          28    76  1881      0.112
28   415    27    28     1      1805 female 3      lower yes          28    77  1882     -0.068
29   463     0     1     0      1836 male   3      upper yes           3    50  1886      0.129
30   463     1     2     0      1836 male   3      upper yes           3    51  1887     -0.378
31   463     2     3     0      1836 male   3      upper yes           3    52  1888     -0.08 
32   179     0     1     0      1784 male   2      lower yes           5    50  1834     -0.177
33   179     1     2     0      1784 male   2      lower yes           5    51  1835     -0.085
34   179     2     3     0      1784 male   2      lower yes           5    52  1836     -0.007
35   179     3     4     0      1784 male   2      lower yes           5    53  1837      0.104
36   179     4     5     1      1784 male   2      lower yes           5    54  1838      0.118
37   526     0     1     0      1777 male   4      lower yes          19    50  1827     -0.037
38   526     1     2     0      1777 male   4      lower yes          19    51  1828     -0.261
39   526     2     3     0      1777 male   4      lower yes          19    52  1829      0.044
40   526     3     4     0      1777 male   4      lower yes          19    53  1830      0.126
41   526     4     5     0      1777 male   4      lower yes          19    54  1831      0.423
42   526     5     6     0      1777 male   4      lower yes          19    55  1832     -0.019
43   526     6     7     0      1777 male   4      lower yes          19    56  1833     -0.156
44   526     7     8     0      1777 male   4      lower yes          19    57  1834     -0.177
45   526     8     9     0      1777 male   4      lower yes          19    58  1835     -0.085
46   526     9    10     0      1777 male   4      lower yes          19    59  1836     -0.007
47   526    10    11     0      1777 male   4      lower yes          19    60  1837      0.104
48   526    11    12     0      1777 male   4      lower yes          19    61  1838      0.118
49   526    12    13     0      1777 male   4      lower yes          19    62  1839     -0.197
50   526    13    14     0      1777 male   4      lower yes          19    63  1840     -0.069
```

A very nice tutoral on how to prepare your data can be found here: https://www.rensvandeschoot.com/tutorials/discrete-time-survival/
Specifically, the R package '''discSurv''' is a great package to prepare your dataset.

To run the model, you first need to source the background functions.

```
source("~/R/G2G/G2G_background_functions.r")
```

For the data example above, you just need to run this command:

```
G2G_varying_MLE(Surv(exit,event) ~ sex + immigrant + foodprices, data=Scania_PersonPeriod_Train, subject="id") 
```

The left-hand-side of formula follows the standard R convention for survival analysis, ```Surv(period, occur)```.  Remember that the second variable ```occur``` should be 0 unless the event of interest _occurs_ at that period, which then takes the value of 1.  The right-hand-side of the formula consists of the time-varying and also time-invarying covariates. The ```data=``` indicates the data frame.  Finally, the ```subject=``` indicates the column representing the _person_, and note that it should be a string.  You should also standardize or at least downscale your covariates to avoid numerical overflow issue.

In this example, the **exit** variable showing 1, 2, 3,..., etc., represents the period.  It is important to have all the periods for each subject until the end of observation for that subject.

The estimation algorithm for ```G2G_varying_MLE``` is just L-BFGS-B using ```optim``` from R.  The output is the usual information from it, plus the standard errors for the parameters and the Wald confidence limits.

Please contact me (kaloklee@gmail.com) if you have any questions or comments.

Some references for the (shifted) Beta-Geometric Distribution:

Fader, Peter S. and Bruce G.S. Hardie (2007), "How to Project Customer Retention," Journal of Interactive Marketing, 31 (1), 76-90.

Lee, Ka Lok, Peter S. Fader, and Bruce G.S. Hardie (2007), “How to Project Patient Persistency,” Foresight: The International Journal of Applied Forecasting, Issue 8, 31-35.
