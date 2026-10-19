# 5a — Regla de Simpson

El programa calcula la integral de la distribución t de Student usando la regla de Simpson y una tolerancia de `0.00001`.

Desde la carpeta `5a`:

```bash
mvn test
mvn compile
java -cp target/classes App
```

Al ejecutarse solicita, en este orden, el valor inicial, el valor final y los grados de libertad (DOF). Muestra el resultado en consola y lo guarda en `output/output5.txt`.
