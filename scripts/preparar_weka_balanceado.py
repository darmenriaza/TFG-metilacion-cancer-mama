import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path("datos")

METHYLATION_FILE = DATA_DIR / "TCGA-BRCA.methylation450.tsv.gz"

# Número de CpGs que se quiere conservar
N_CPGS = 1000

# Semilla para que el muestreo sea reproducible
RANDOM_STATE = 42

# Lectura por bloques para no cargar todo de golpe
CHUNK_SIZE = 2000

OUTPUT_CSV = f"weka_brca_balanceado_{N_CPGS}_cpgs.csv"
OUTPUT_ARFF = f"weka_brca_balanceado_{N_CPGS}_cpgs.arff"
OUTPUT_SAMPLES = f"sample_ids_balanceado_{N_CPGS}_cpgs.csv"


def get_sample_type(sample_id):
    """
    Extrae el código de tipo de muestra del barcode TCGA.
    Ejemplo:
    TCGA-AQ-A0Y5-01A -> 01 = tumor primario
    TCGA-BH-A0B3-11A -> 11 = tejido normal sólido
    """
    parts = sample_id.split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return None



# Leer cabecera y seleccionar muestras

print("Leyendo cabecera del archivo de metilación...")

header = pd.read_csv(METHYLATION_FILE, sep="\t", nrows=0).columns.tolist()

id_col = header[0]
samples = header[1:]

tumor_samples = [s for s in samples if get_sample_type(s) == "01"]
normal_samples = [s for s in samples if get_sample_type(s) == "11"]

print("Muestras tumorales disponibles:", len(tumor_samples))
print("Muestras normales disponibles:", len(normal_samples))

# Submuestreo: seleccionar el mismo número de tumores que normales
np.random.seed(RANDOM_STATE)
tumor_samples_balanced = np.random.choice(
    tumor_samples,
    size=len(normal_samples),
    replace=False
).tolist()

selected_samples = tumor_samples_balanced + normal_samples

labels = pd.Series(
    {
        s: "Tumor" if get_sample_type(s) == "01" else "Normal"
        for s in selected_samples
    },
    name="class"
)

print()
print("Dataset balanceado:")
print(labels.value_counts())
print("Total de muestras:", len(selected_samples))
print()


# Selección de CpGs por variabilidad


print("Calculando la variabilidad de las CpGs...")

usecols = [id_col] + selected_samples
variances = []

reader = pd.read_csv(
    METHYLATION_FILE,
    sep="\t",
    usecols=usecols,
    chunksize=CHUNK_SIZE
)

for i, chunk in enumerate(reader, start=1):
    chunk = chunk.set_index(id_col)

    # Convertir valores a numéricos
    chunk = chunk.apply(pd.to_numeric, errors="coerce")

    # Eliminar CpGs con demasiados valores perdidos
    min_valid = int(0.8 * len(selected_samples))
    chunk = chunk.dropna(thresh=min_valid)

    # Calcular varianza de cada CpG
    chunk_var = chunk.var(axis=1, skipna=True)
    variances.append(chunk_var)

    if i % 20 == 0:
        print(f"Bloques procesados: {i}")

all_variances = pd.concat(variances)

top_cpgs = all_variances.nlargest(N_CPGS).index.tolist()

print()
print(f"CpGs seleccionadas: {len(top_cpgs)}")
print("Primeras CpGs seleccionadas:")
print(top_cpgs[:10])
print()


# Extraer las CpGs seleccionadas


print("Extrayendo las CpGs seleccionadas...")

selected_chunks = []
top_cpgs_set = set(top_cpgs)

reader = pd.read_csv(
    METHYLATION_FILE,
    sep="\t",
    usecols=usecols,
    chunksize=CHUNK_SIZE
)

for i, chunk in enumerate(reader, start=1):
    chunk = chunk.set_index(id_col)

    # Quedarse solo con las CpGs seleccionadas
    chunk = chunk.loc[chunk.index.intersection(top_cpgs_set)]

    if not chunk.empty:
        chunk = chunk.apply(pd.to_numeric, errors="coerce")
        selected_chunks.append(chunk)

    if i % 20 == 0:
        print(f"Bloques procesados: {i}")

meth_top = pd.concat(selected_chunks)

# Reordenar según la lista de CpGs más variables
meth_top = meth_top.loc[top_cpgs]

print("Matriz CpG x muestras:", meth_top.shape)


# Preparar matriz para Weka

print("Preparando matriz para Weka...")

# Transponer: filas = muestras, columnas = CpGs
X = meth_top.T

# Imputar valores perdidos con la mediana de cada CpG
X = X.fillna(X.median(axis=0))

# Añadir clase como última columna
X["class"] = labels.loc[X.index]

# Guardar IDs de muestras aparte para trazabilidad
sample_ids = pd.DataFrame({
    "sample_id": X.index,
    "class": X["class"].values
})
sample_ids.to_csv(OUTPUT_SAMPLES, index=False)

# Quitar el índice para que Weka no use el ID como variable
X_for_weka = X.reset_index(drop=True)

# Guardar CSV
X_for_weka.to_csv(OUTPUT_CSV, index=False)

print(f"Archivo CSV generado: {OUTPUT_CSV}")
print(f"Archivo con IDs de muestras generado: {OUTPUT_SAMPLES}")


# Crear archivo ARFF para Weka


with open(OUTPUT_ARFF, "w", encoding="utf-8") as f:
    f.write("@RELATION tcga_brca_methylation_balanced\n\n")

    for col in X_for_weka.columns:
        if col != "class":
            f.write(f"@ATTRIBUTE {col} NUMERIC\n")

    f.write("@ATTRIBUTE class {Normal,Tumor}\n\n")
    f.write("@DATA\n")

    for _, row in X_for_weka.iterrows():
        values = []
        for col in X_for_weka.columns:
            values.append(str(row[col]))
        f.write(",".join(values) + "\n")

print(f"Archivo ARFF generado: {OUTPUT_ARFF}")
