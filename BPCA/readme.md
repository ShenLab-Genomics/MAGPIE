BPCAfill implementation used in MAGPIE is modified from original repository: [bpcafill](https://github.com/shigeyukioba/bpcafill/tree/master). 
The usage of BPCAfill for MAGPIE is identical to the original one. 

This is the guidance of original repo.

```
Prepare data matrix in which missing entries in matrix are denoted by value 999.0. 
BPCAfill estimates the missing entries using BPCA model and returns complete matrix filled with the estimated values.
>>>x999 = load('sample999.dat');
>>>xfilled = BPCAfill( x999 );
epoch=10, dtau=1.25922
epoch=20, dtau=0.0970274
epoch=30, dtau=0.0349725
epoch=40, dtau=0.0150597
epoch=50, dtau=0.00761689
epoch=60, dtau=0.0034555
epoch=70, dtau=0.000254066
epoch=80, dtau=4.42614e-005    <-- when dtau < 1e-004, it is decided as converged.
>> save FilledExpression.dat -ascii xfilled
Please wait for the calculation to converge. It will take some minutes or hours till finish. 
(Example. The sample999.dat (757x50 entries with 797 missing) costs 10 minutes to fill using Pentium III 600MHz processer. )
```
If interested, please refer to [original repo](https://github.com/shigeyukioba/bpcafill/tree/master) for further usage guidance.