# Matrix_multiplication

* compile
```
make
```

* run
```
make run
// just show how long it take.
```

* generate test data
```
./gen_test 1024 1024 > input.txt
```

* compute
```
./matrix_mul < input.txt
./matrix_Strassen < input.txt
```

## note
`display(A)` will show all elements of A.

* matrix_Strassen 
  * implement of Strassen algorithm
  
* matrix_Strassen_parallel 
  * using omp to paralize
  
* matrix_mul 
  * traditional matrix multiplication with pthread
  
* matrix_opt 
  * cache friendly by using transpose matrix
