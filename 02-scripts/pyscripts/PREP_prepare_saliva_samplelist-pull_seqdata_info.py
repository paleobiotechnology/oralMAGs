import re

import pandas as pd
from amdirt.core import ena

publication_year = re.search(r'[A-Za-z]+([0-9]+)',
                             snakemake.wildcards.study).group(1)

e = ena.ENAPortalAPI()
ena_results = e.query(snakemake.params.accession,
                      fields=["study_accession",
                              "run_accession",
                              "sample_alias",
                              "secondary_sample_accession",
                              "library_layout",
                              "fastq_ftp",
                              "fastq_md5",
                              "tax_id",
                              "instrument_model",
                              "read_count"])
seqdata = pd.DataFrame(ena_results)
if snakemake.wildcards.study == "Clemente2015":
    seqdata = seqdata.query('sample_alias.str.startswith("O")')
else:
    seqdata = seqdata.query('tax_id == "256318"')
seqdata['project_name'] = snakemake.wildcards.study
seqdata['publication_year'] = publication_year

(seqdata
    .iloc[:, [10, 11] + list(range(0, 10))]
    .to_csv(snakemake.output[0], sep="\t", index=False)
)
