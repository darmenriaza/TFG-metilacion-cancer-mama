import random

input_file = "weka_brca_balanceado_1000_cpgs.arff"
output_file = "weka_brca_balanceado_1000_cpgs_clases_mezcladas.arff"

random.seed(42)

with open(input_file, "r", encoding="utf-8") as f:
    lines = f.readlines()

header = []
data = []
in_data = False

for line in lines:
    stripped = line.strip()

    if stripped.lower() == "@data":
        header.append(line)
        in_data = True
        continue

    if not in_data:
        header.append(line)
    else:
        if stripped == "" or stripped.startswith("%"):
            data.append(line)
        else:
            data.append(line)

# Separar las filas reales de datos
data_rows = [line.strip() for line in data if line.strip() and not line.strip().startswith("%")]

# Extraer clases
features = []
classes = []

for row in data_rows:
    parts = row.split(",")
    features.append(parts[:-1])
    classes.append(parts[-1])

# Mezclar solo las etiquetas
shuffled_classes = classes.copy()
random.shuffle(shuffled_classes)

# Crear nuevas filas
new_data_rows = []
for feat, new_class in zip(features, shuffled_classes):
    new_data_rows.append(",".join(feat + [new_class]) + "\n")

with open(output_file, "w", encoding="utf-8") as f:
    f.writelines(header)
    f.writelines(new_data_rows)

print("Archivo generado:", output_file)
print("Clases originales:", {c: classes.count(c) for c in set(classes)})
print("Clases mezcladas:", {c: shuffled_classes.count(c) for c in set(shuffled_classes)})