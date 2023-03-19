Optims:
1. usar Static Vectors
2. usar una tabla de las entradas `J`, `K` para no recalcular
3. usar `FillArrays` para los `\alpha`
4. NO USAR UN MALDITO TRY CATCH
5. `unique!` could be run on a `Set`, no need to go over the array twice!