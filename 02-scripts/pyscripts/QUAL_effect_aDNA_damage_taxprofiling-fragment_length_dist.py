import json
import os

import pandas as pd

fns = snakemake.input.modern + snakemake.input.ancient

results = []
for fn in fns:
    sample = os.path.basename(fn).replace(".fastp.json", "")
    with open(fn, "rt") as jsonfile:
        fastp = json.load(jsonfile)
    histogram = pd.DataFrame.from_dict({'counts': fastp['insert_size']['histogram']}).transpose()
    histogram.index = [sample]
    histogram['unknown'] = fastp['insert_size']['unknown']
    results.append(histogram)

df = pd.concat(results)
df.columns = list(range(1, 513)) + ["unknown"]
df.to_csv(snakemake.output[0], sep="\t", index=True)
