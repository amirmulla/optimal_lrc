## Erasure Test

```
GF:  13
Code dim n:  12
Message dim k:  8
Locality r:  3
Global Minimum distance d:  3
Local Minimum distance d:  2
[12, 8, 3] Tamo-Berg Code over GF(13)
Erasure Vector:  (1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0)
Recieved Word :  (0, 5, 7, 4, 9, 2, 9, 0, 0, 11, 5, 7)
Codeword      :  (7, 5, 7, 4, 9, 2, 9, 0, 11, 11, 5, 7)
Corr Codeword :  (7, 5, 7, 4, 9, 2, 9, 0, 11, 11, 5, 7)
Correction Successfull:  True
```

## Error Test

```
GF:  16
Code dim n:  16
Message dim k:  5
Locality r:  2
Local Minimum distance d:  3
[16, 5, 2] Tamo-Berg Code over GF(16)
Error              :  (0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 15, 0)
Recieved Word      :  (12, 4, 9, 1, 8, 5, 0, 13, 5, 6, 5, 5, 4, 1, 6, 12)
Original Codeword  :  (12, 4, 9, 1, 8, 5, 0, 13, 5, 5, 5, 5, 4, 1, 9, 12)
Corrected Codeword :  (12, 4, 9, 1, 8, 5, 0, 13, 5, 5, 5, 5, 4, 1, 9, 12)
Correction Successfull:  True
```

## TwoSets Error Test

```
GF:  64
Code dim n:  63
Message dim k:  35
Locality r1:  5
Locality r2:  7
Local Minimum distance d1:  3
Local Minimum distance d2:  3
[63, 35, [5, 7]] Tamo-Berg Code over GF(64)
Error              :  (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 29, 0, 40, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 38, 0, 0, 0, 0, 21, 0, 8)
Recieved Word      :  (21, 44, 18, 44, 53, 53, 19, 28, 25, 41, 44, 34, 43, 3, 35, 13, 9, 31, 2, 23, 39, 39, 42, 47, 30, 37, 49, 32, 26, 62, 9, 38, 59, 51, 22, 49, 39, 42, 19, 32, 49, 20, 2, 20, 52, 2, 8, 42, 36, 50, 43, 55, 31, 23, 3, 62, 60, 54, 56, 55, 10, 33, 1)
Original Codeword  :  (21, 44, 18, 44, 53, 53, 19, 28, 25, 41, 44, 34, 43, 3, 35, 16, 9, 55, 2, 23, 39, 39, 42, 47, 30, 37, 49, 32, 26, 62, 9, 38, 59, 51, 22, 49, 39, 42, 19, 32, 49, 20, 2, 20, 52, 2, 8, 42, 36, 50, 43, 55, 31, 23, 3, 24, 60, 54, 56, 55, 31, 33, 9)
Corrected Codeword :  ((21, 44, 18, 44, 53, 53, 19, 28, 25, 41, 44, 34, 43, 3, 35, 16, 9, 55, 2, 23, 39, 39, 42, 47, 30, 37, 49, 32, 26, 62, 9, 38, 59, 51, 22, 49, 39, 42, 19, 32, 49, 20, 2, 20, 52, 2, 8, 42, 36, 50, 43, 55, 31, 23, 3, 24, 60, 54, 56, 55, 31, 33, 9), 2)
Correction Successfull:  True
```