import pandas as pd
import numpy as np
from pathlib import Path

DATA_DIR = Path("datos")

METHYLATION_FILE = DATA_DIR / "TCGA-BRCA.methylation450.tsv.gz"
EVOLUTION_FILE = Path("pacientes_evolucion_5_anios.csv")

N_CPGS = 1000
CHUNK_SIZE = 2000

OUTPUT_CSV = f"weka_brca_evolucion_5_anios_{N_CPGS}_cpgs.csv"
OUTPUT_ARFF = f"weka_brca_evolucion_5_anios_{N_CPGS}_cpgs.arff"
OUTPUT_SAMPLES = f"sample_ids_evolucion_5_anios_{N_CPGS}_cpgs.csv"


print("Leyendo pacientes con evolución clínica...")

evol = pd.read_csv(EVOLUTION_FILE)

evol = evol.dropna(subset=["evolucion_5_anios"]).copy()

selected_samples = evol["sample_id"].astype(str).tolist()

labels = pd.Series(
    evol["evolucion_5_anios"].values,
    index=evol["sample_id"].astype(str),
    name="class"
)

print("Pacientes útiles:", len(selected_samples))
print(labels.value_counts())
print()


print("Leyendo cabecera del archivo de metilación...")

header = pd.read_csv(METHYLATION_FILE, sep="\t", nrows=0).columns.tolist()

id_col = header[0]
available_samples = set(header[1:])

selected_samples = [s for s in selected_samples if s in available_samples]
labels = labels.loc[selected_samples]

print("Muestras encontradas en metilación:", len(selected_samples))
print(labels.value_counts())
print()



# Selección de CpGs por varianza

print("Calculando varianza de las CpGs...")

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

    chunk = chunk.apply(pd.to_numeric, errors="coerce")

    # Conservar CpGs con al menos 80 % de valores válidos
    min_valid = int(0.8 * len(selected_samples))
    chunk = chunk.dropna(thresh=min_valid)

    # Calcular varianza ignorando ausentes
    chunk_var = chunk.var(axis=1, skipna=True)

    variances.append(chunk_var)

    if i % 20 == 0:
        print(f"Bloques procesados: {i}")

all_variances = pd.concat(variances)

top_cpgs = all_variances.nlargest(N_CPGS).index.tolist()

print()
print("CpGs seleccionadas:", len(top_cpgs))
print("Primeras CpGs seleccionadas:")
print(top_cpgs[:10])
print()



# Extraer CpGs seleccionadas

print("Extrayendo CpGs seleccionadas...")

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

    chunk = chunk.loc[chunk.index.intersection(top_cpgs_set)]

    if not chunk.empty:
        chunk = chunk.apply(pd.to_numeric, errors="coerce")
        selected_chunks.append(chunk)

    if i % 20 == 0:
        print(f"Bloques procesados: {i}")

meth_top = pd.concat(selected_chunks)

meth_top = meth_top.loc[top_cpgs]

print("Matriz CpG x muestras:", meth_top.shape)



# Preparar matriz para Weka

print("Preparando matriz para Weka...")

X = meth_top.T

# Imputar valores ausentes restantes con la mediana de cada CpG
X = X.fillna(X.median(axis=0))

X["class"] = labels.loc[X.index]

sample_ids = pd.DataFrame({
    "sample_id": X.index,
    "class": X["class"].values
})
sample_ids.to_csv(OUTPUT_SAMPLES, index=False)

X_for_weka = X.reset_index(drop=True)

X_for_weka.to_csv(OUTPUT_CSV, index=False)

print("Archivo CSV generado:", OUTPUT_CSV)
print("Archivo con IDs generado:", OUTPUT_SAMPLES)



# Crear ARFF


print("Generando archivo ARFF...")

with open(OUTPUT_ARFF, "w", encoding="utf-8") as f:
    f.write("@RELATION tcga_brca_evolucion_5_anios\n\n")

    for col in X_for_weka.columns:
        if col != "class":
            f.write(f"@ATTRIBUTE {col} NUMERIC\n")

    f.write("@ATTRIBUTE class {Buena_evolucion,Mala_evolucion}\n\n")
    f.write("@DATA\n")

    for _, row in X_for_weka.iterrows():
        values = []
        for col in X_for_weka.columns:
            values.append(str(row[col]))
        f.write(",".join(values) + "\n")

print("Archivo ARFF generado:", OUTPUT_ARFF)
