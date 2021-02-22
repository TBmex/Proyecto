## Baps Valencia
- El programa utilizado es **Fastbaps** el cual parte de un alineamiento previamente obtenido por el pipeline  
- El alineamiento usado en este caso sera un alineamiento reducido de 1177 de linaje cuatro de genomas de pacientes de TB  
- Se reducira siguiendo la siguiente condición: Un solo genoma por cluster español

### Obteniendo alineamiento de solo ESP
input: **run_alignment_no_resis_1177.fas**  
Utilizares el **scrip seqtk** para darle una lista de aislados reducida

~~~

Ejemplo
Extract sequences with names in file name.lst, one sequence name per line:

  > seqtk subseq in.fq name.lst > out.fq

Aplicacion

  > seqtk subseq run_alignment_no_resis_1177.fas ID_ESP.lst > ESP.fq

~~~
Creamos la lista ID_ESP.lst filtrando con Exel, se requiere codigo para agilizar esta parte.
~~~
Crear codigo para obtener una muestra por cluster en tablas de R.
~~~

Conclusión se obtuvo de output el archivo: **ESP.fq**

### Generando tabla de Genotipos con Fastbaps
input: **ESP.fq**  
Utilizamos el scrip llamado **"Proyecto"**
~~~
#Libraries
library(fastbaps)
library(ape)

#Loading data OK
fasta.file.name <- "ESP.fq"
sparse.data <- import_fasta_sparse_nt(fasta.file.name)

#Detalle OK
sparse.data <- optimise_prior(sparse.data, type = "baps")

#Running fastbaps "baps.hc" es el archivo a obtener.
baps.hc <- fast_baps(sparse.data)

# Obtener tablas

#Bayesian hierarchical partition
ESP_baps <- multi_res_baps(sparse.data, levels = 4)
~~~
Resultados

~~~
> head(ESP_baps)
  Isolates Level 1 Level 2 Level 3 Level 4
1      G01       1       2       4      13
2      G02       1       2       4      13
3      G03       1       2       4      14
4    G1000       1       2       5      16
5    G1002       2       6      16      43
6    G1009       4      14      40     110
~~~
Conclusión se obtuvo de output la tabla: **ESP_baps**  
Nota: tambien hay un valor de bootstrap

### Generacion de arbol filogenetico (ITOL)
input: **ESP.fq**

Mejor un arbol nuevo o sobreponerlo?
Corriendo el arbol...
~~~
/home/carlos/Epi_val_baps

nohup nice -n 5 iqtree -s ESP.fq -m GTR -nt 20 -o G177m -vv -bb 1000 &
~~~


### Tablas de frecuencias y asociaciones
