import pandas as pd
import numpy as np
from pathlib import Path
from collections import Counter

DATA_DIR = Path("datos")

METHYLATION_FILE = DATA_DIR / "TCGA-BRCA.methylation450.tsv.gz"
CLINICAL_FILE = DATA_DIR / "TCGA-BRCA.clinical.tsv.gz"


def get_sample_type(sample_id):
    parts = str(sample_id).split("-")
    if len(parts) >= 4:
        return parts[3][:2]
    return None


def get_patient_id(barcode):
    return str(barcode)[:12]


def clean_numeric(series):
    return pd.to_numeric(series, errors="coerce")


def parse_event_from_vital(series):
    def parse_value(x):
        if pd.isna(x):
            return np.nan

        s = str(x).strip().lower()

        if s in ["1", "1.0"]:
            return 1
        if s in ["0", "0.0"]:
            return 0

        if "dead" in s or "deceased" in s:
            return 1
        if "alive" in s or "living" in s:
            return 0

        return np.nan

    return series.apply(parse_value)


def find_column(columns, exact_names):
    lower_map = {c.lower(): c for c in columns}
    for name in exact_names:
        if name.lower() in lower_map:
            return lower_map[name.lower()]
    return None



# Leer muestras tumorales del archivo de metilación

print("Leyendo cabecera del archivo de metilación...")

header = pd.read_csv(METHYLATION_FILE, sep="\t", nrows=0).columns.tolist()
samples = header[1:]

tumor_samples = [s for s in samples if get_sample_type(s) == "01"]

tumor_df = pd.DataFrame({
    "sample_id": tumor_samples
})
tumor_df["patient_id"] = tumor_df["sample_id"].apply(get_patient_id)

print("\nMUESTRAS DE METILACIÓN")
print("Muestras tumorales primarias 01:", len(tumor_df))
print("Pacientes únicos con muestra tumoral 01:", tumor_df["patient_id"].nunique())


# Leer archivo clínico

print("\nLeyendo archivo clínico...")

clinical = pd.read_csv(CLINICAL_FILE, sep="\t", low_memory=False)

print("\nARCHIVO CLÍNICO")
print("Dimensiones:", clinical.shape)
print("\nColumnas del archivo clínico:")
for col in clinical.columns:
    print("-", col)



# Mostrar columnas candidatas relacionadas con evolución

keywords = [
    "vital", "death", "dead", "follow", "surviv", "survival",
    "days", "os", "dss", "dfi", "pfi", "progress", "recur",
    "relapse", "status"
]

candidate_cols = [
    c for c in clinical.columns
    if any(k in c.lower() for k in keywords)
]

print("\nCOLUMNAS CANDIDATAS RELACIONADAS CON EVOLUCIÓN/SUPERVIVENCIA")
for col in candidate_cols:
    print("-", col)



# Detectar columna identificadora de paciente o muestra


possible_id_cols = [
    "sample",
    "sample_id",
    "submitter_id",
    "bcr_patient_barcode",
    "patient",
    "patient_id",
    "case_submitter_id"
]

id_col = find_column(clinical.columns, possible_id_cols)

if id_col is None:
    # Buscar automáticamente una columna con barcodes TCGA
    for col in clinical.columns:
        values = clinical[col].dropna().astype(str)
        if len(values) > 0:
            prop_tcga = values.str.startswith("TCGA-").mean()
            if prop_tcga > 0.5:
                id_col = col
                break

if id_col is None:
    print("\nNo se ha podido detectar una columna identificadora TCGA.")
    print("Revisa las columnas impresas arriba y dime cuál parece contener los IDs de paciente o muestra.")
    raise SystemExit

print("\nColumna identificadora detectada:", id_col)

clinical["patient_id"] = clinical[id_col].apply(get_patient_id)

clinical_patients = clinical.drop_duplicates(subset="patient_id").copy()

merged = tumor_df.merge(
    clinical_patients,
    on="patient_id",
    how="left",
    indicator=True
)

print("\nUNIÓN METILACIÓN + CLÍNICA")
print("Pacientes tumorales en metilación:", tumor_df["patient_id"].nunique())
print("Pacientes enlazados con clínica:", (merged["_merge"] == "both").sum())
print("Pacientes sin clínica encontrada:", (merged["_merge"] == "left_only").sum())



# Intentar construir variable de supervivencia/evolución

# si existen columnas tipo OS y OS.time
os_col = find_column(merged.columns, ["OS", "overall_survival", "overall_survival_event"])
os_time_col = find_column(merged.columns, ["OS.time", "OS_time", "overall_survival_time"])

# si existen columnas clínicas básicas
vital_col = find_column(merged.columns, [
    "vital_status",
    "vital_status.demographic",
    "vital_status.demographic.vital_status"
])

days_death_col = find_column(merged.columns, [
    "days_to_death",
    "days_to_death.demographic",
    "days_to_death.diagnoses"
])

days_follow_col = find_column(merged.columns, [
    "days_to_last_follow_up",
    "days_to_last_followup",
    "days_to_last_follow_up.diagnoses",
    "days_to_last_followup.diagnoses",
    "days_to_last_known_alive",
    "days_to_last_known_alive.demographic"
])

analysis = merged.drop_duplicates(subset="patient_id").copy()

event = None
time = None

if os_col is not None and os_time_col is not None:
    print("\nSe usarán columnas OS ya procesadas:")
    print("Evento:", os_col)
    print("Tiempo:", os_time_col)

    event = parse_event_from_vital(analysis[os_col])
    time = clean_numeric(analysis[os_time_col])

elif vital_col is not None and (days_death_col is not None or days_follow_col is not None):
    print("\nSe construirán variables de supervivencia a partir de columnas clínicas:")
    print("Estado vital:", vital_col)
    print("Días hasta fallecimiento:", days_death_col)
    print("Días hasta último seguimiento:", days_follow_col)

    event = parse_event_from_vital(analysis[vital_col])

    days_death = clean_numeric(analysis[days_death_col]) if days_death_col is not None else pd.Series(np.nan, index=analysis.index)
    days_follow = clean_numeric(analysis[days_follow_col]) if days_follow_col is not None else pd.Series(np.nan, index=analysis.index)

    time = days_death.copy()
    time = time.fillna(days_follow)

else:
    print("\nNo se han encontrado automáticamente columnas suficientes para construir supervivencia.")
    print("Pásame la lista de columnas candidatas que aparece arriba.")
    raise SystemExit


analysis["evento_muerte"] = event
analysis["tiempo_seguimiento_dias"] = time

print("\nRESUMEN DE SUPERVIVENCIA")
print("Pacientes con evento conocido:", analysis["evento_muerte"].notna().sum())
print("Pacientes con tiempo de seguimiento conocido:", analysis["tiempo_seguimiento_dias"].notna().sum())

print("\nEvento de muerte:")
print(analysis["evento_muerte"].value_counts(dropna=False).rename(index={0: "Vivo/censurado", 1: "Fallecido"}))


# Variable supervivencia a 5 años


limite_5_anios = 5 * 365.25

def classify_5_year(row):
    event = row["evento_muerte"]
    time = row["tiempo_seguimiento_dias"]

    if pd.isna(event) or pd.isna(time):
        return np.nan

    # Mala evolución: fallece antes de 5 años
    if event == 1 and time < limite_5_anios:
        return "Mala_evolucion"

    # Buena evolución: alcanza al menos 5 años de seguimiento o supervivencia
    if time >= limite_5_anios:
        return "Buena_evolucion"

    # Vivo/censurado con menos de 5 años: no sabemos si alcanzará los 5 años
    return np.nan


analysis["evolucion_5_anios"] = analysis.apply(classify_5_year, axis=1)

print("\nCLASIFICACIÓN EXPLORATORIA: SUPERVIVENCIA A 5 AÑOS")
print(analysis["evolucion_5_anios"].value_counts(dropna=False))

usable = analysis.dropna(subset=["evolucion_5_anios"])

print("\nPacientes útiles para un modelo binario de evolución a 5 años:", len(usable))
print(usable["evolucion_5_anios"].value_counts())

# Guardar tabla para revisarla
output_file = "pacientes_evolucion_5_anios.csv"
analysis[[
    "sample_id",
    "patient_id",
    "evento_muerte",
    "tiempo_seguimiento_dias",
    "evolucion_5_anios"
]].to_csv(output_file, index=False)

print("\nArchivo generado:", output_file)
