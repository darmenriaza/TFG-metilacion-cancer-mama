# Datos utilizados en Weka

Esta carpeta contiene los conjuntos de datos procesados utilizados en los experimentos realizados con Weka.

## weka_brca_balanceado_1000_cpgs.arff

Conjunto de datos utilizado para la clasificación entre muestras tumorales y muestras normales de cáncer de mama.  
Está formado por 194 muestras balanceadas, 97 tumorales y 97 normales, y contiene las 1000 posiciones CpG seleccionadas durante el preprocesamiento.

## weka_brca_balanceado_1000_cpgs_clases_mezcladas.arff

Conjunto de datos utilizado como comprobación adicional.  
Mantiene los mismos valores de metilación que el conjunto balanceado, pero con las etiquetas de clase mezcladas aleatoriamente. Se utilizó para comprobar si los modelos seguían obteniendo buenos resultados al romper la relación real entre los perfiles de metilación y la clase.

## weka_brca_evolucion_5_anios_1000_cpgs.arff

Conjunto de datos utilizado para el análisis exploratorio de evolución clínica a 5 años.  
Incluye únicamente pacientes con evolución conocida y contiene 1000 posiciones CpG como variables de entrada, junto con una clase final que indica buena o mala evolución.

## Nota

Los archivos originales de TCGA-BRCA no se incluyen en este repositorio debido a su tamaño. Los conjuntos incluidos corresponden a versiones procesadas utilizadas en los experimentos del TFG.
