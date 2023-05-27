# Gram–Schmidt process
Реализация метода ортогонализации [Грама-Шмидта](https://ru.wikipedia.org/wiki/Процесс_Грама_―_Шмидта) с внедрением КЭШ оптимизации, получаемой за счёт использования блочности алгоритма.
Алгоритм был написан, в нескольких версиях:
- Only Fortran
- Only C++
- BLAS Fortran
- BLAS C++
- LAPACK Fortran
- LAPACK C++

Результаты численных экспериментов приведены в [отчёте](https://github.com/DrEternity/Gram_Schmidt_process/blob/main/presentation.pdf)

![Псевдокод алгоритма:](https://github.com/DrEternity/Gram_Schmidt_process/blob/main/pseudocode.jpg?raw=true)
