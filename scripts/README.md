# Scripts utilizados en el TFG

Esta carpeta contiene los scripts principales utilizados para preparar los datos empleados en los experimentos realizados en Weka.

## preparar_weka_balanceado.py

Este script genera el conjunto de datos utilizado en el experimento de clasificación entre muestras tumorales y muestras normales.  
A partir del archivo de metilación de TCGA-BRCA, selecciona las muestras tumorales primarias y las muestras normales sólidas, balancea las clases, selecciona las 1000 posiciones CpG con mayor varianza, imputa los valores ausentes mediante la mediana de cada CpG y genera los archivos finales en formato CSV y ARFF para Weka.

## comprobar_evolucion_clinica.py

Este script se utilizó para construir la variable de evolución clínica a 5 años.  
Para ello, combina las muestras tumorales del archivo de metilación con la información clínica disponible, identifica el estado vital y el tiempo de seguimiento de cada paciente, y clasifica a los pacientes en dos grupos: buena evolución y mala evolución. El resultado se guarda en el archivo `pacientes_evolucion_5_anios.csv`.

## crear_dataset_evolución.py

Este script genera el conjunto de datos utilizado en el experimento de predicción exploratoria de evolución clínica a 5 años.  
Utiliza el archivo `pacientes_evolucion_5_anios.csv`, selecciona los pacientes con evolución conocida, extrae sus perfiles de metilación, selecciona las 1000 posiciones CpG con mayor varianza, imputa los valores ausentes mediante la mediana y crea los archivos CSV y ARFF para Weka.

## mezclar_clases.py

Este script se utilizó como comprobación adicional.  
A partir del conjunto balanceado de clasificación tumoral, mantiene los valores de metilación originales pero mezcla aleatoriamente las etiquetas de clase. De esta forma se genera un conjunto en el que se rompe la relación real entre los perfiles de metilación y la clase, permitiendo comprobar si los modelos siguen obteniendo buenos resultados cuando las etiquetas no contienen información real.
