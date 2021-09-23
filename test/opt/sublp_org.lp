minimize
obj: 1 x1 + 1e3 x3 + 1e3 x4 + 1e3 x5 + 1e3 x6 + 1e3 x7 + 1e3 x8 + 1e3 x9 + 1e3 x10
subject to
c3: -1 x7 + 1 x3 + 1 x1 >= -2
c4: 1 x3 >= 0
c5: 1 x7 >= 0
c6: 1 x4 >= 0
c7: 1 x8 >= 0
c8: 1 x5 >= 0
c9: 1 x9 >= 0
c10: -1 x10 + 1 x6 >= 0
c11: 1 x6 >= 0
c12: 1 x10 >= 0
c1: -1 x8 + 1 x4 - 1 x1 = 2
c2: -1 x9 + 1 x5 = 1
Bounds
x1 <= 1e3
x2 <= 1e3
x1 >= -1e3
x2 >= -1e3
x3 free
x4 free
x5 free
x6 free
x7 free
x8 free
x9 free
x10 free
End
