# Repetición de mini prueba

En el código modifiqué la función *get_all_metadata_from_GSE()* para obtener los archivos suplementarios de los GSMs antes de mergear su DataFrame con el del GSE y los GPls. QUiero ver que no modifique el resultado antes de intentar hacer la prueba masiva.

```Bash
nohup python MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py -t TxIDs.txt -o 23082023Test.tsv -e diegocar@lcg.unam.mx -a 34677fdcfd2f0659a7f9ee05ab6e44704f09 &
```

Nota: hice ligas simbólicas a los archivos 

No hubo errores

# Repetición de la prueba grande

```Bash
nohup python MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py -t abasy.txt -o 27082023Test.tsv -e diegocar@lcg.unam.mx -a 34677fdcfd2f0659a7f9ee05ab6e44704f09 &
```

Obtuve este error
```Python
Traceback (most recent call last):
  File "MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py", line 225, in <module>
    metadata_async = list(executor.map(get_GEO_metadata_from_txid, GSE_IDs))
  File "/export/apps/bioconda/lib/python3.8/concurrent/futures/_base.py", line 611, in result_iterator
    yield fs.pop().result()
  File "/export/apps/bioconda/lib/python3.8/concurrent/futures/_base.py", line 432, in result
    return self.__get_result()
  File "/export/apps/bioconda/lib/python3.8/concurrent/futures/_base.py", line 388, in __get_result
    raise self._exception
  File "/export/apps/bioconda/lib/python3.8/concurrent/futures/thread.py", line 57, in run
    result = self.fn(*self.args, **self.kwargs)
  File "MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py", line 189, in get_GEO_metadata_from_txid
    GSE_meta_list = [get_all_metadata_from_GSE(gse, taxonID) for gse in GSE_list]
  File "MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py", line 189, in <listcomp>
    GSE_meta_list = [get_all_metadata_from_GSE(gse, taxonID) for gse in GSE_list]
  File "MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py", line 145, in get_all_metadata_from_GSE
    gse_metadata = get_GEO_meatadata(gse, columns, {'taxonID':txID})
  File "MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py", line 72, in get_GEO_meatadata
    geotype = geo.name[:3]
TypeError: 'NoneType' object is not subscriptable
```

# Segundo intento prueba grande
Resolví el problema agregando un try except

```Bash
nohup python MicroarrayGEOPipeline-ConsultaPrueba_v0.6.py -t abasy.txt -o 27082023Test.tsv -e diegocar@lcg.unam.mx -a 34677fdcfd2f0659a7f9ee05ab6e44704f09 &
```
No hubo errores, commit estable: commit f7577fa9fa1e18d78478c9f815abe435c39cd276

Si el mal de este mundo es no entender los pensamientos del otro, ayúale al prójimo diciéndole lo que piensas.
Diego 2023

# Filtered prueba
Después de agregar función de apply_filters para eliminar Super Series, GPL != Expression profiling by array y txIDs espúreos, hice una corrida de prueba

