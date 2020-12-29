NAME          TESTPROB
ROWS
 N  COST
 E  INF
 L  CON1
 G  CON2
COLUMNS
    X      INF                  3   CON1                 1
    X      CON2                 1   COST                 1
    Y      INF                  8   CON1                 1
    Y      CON2                 1   COST                 1
RHS
    RHS1      CON1                 37   CON1                5
    LHS1      CON2                 0
BOUNDS
 UI BND1      X                 5
 LI BND1      X                 0
 LI BND1      Y                 0
 UI BND1      Y                 5
ENDATA