Optims:
1. usar Static Vectors
2. usar una tabla de las entradas `J`, `K` para no recalcular
3. usar `FillArrays` para los `\alpha`
4. NO USAR UN MALDITO TRY CATCH
5. `unique!` could be run on a `Set`, no need to go over the array twice! Filter before you push anyways
6. Aritmetica de `four_regions` puede precalcular
	- `Ej, Ek, EjOrt, EkOrt, AreaJK, invAreaJK` 
	- (A/_AreaJK_) * dot(_C_, _B_) - (B / _AreaJK_) * dot(_D_, _E_) puede ser 
	- (_invAreaJK_) * (A * dot(_C_, _B_) - B * dot(_D_, _E_)
7. `norm` esta muy caliente - no puedes ser `(e_ij - Site_ij) ^2` < RadiusCluster^2`? no necesita el `sqrt`!