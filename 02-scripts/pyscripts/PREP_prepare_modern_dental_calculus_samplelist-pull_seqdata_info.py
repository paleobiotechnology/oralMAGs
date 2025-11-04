import pandas as pd
from amdirt.core import ena

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
seqdata = (pd.DataFrame(ena_results)
    .query('sample_alias.str.startswith(@snakemake.wildcards.study)')
)
seqdata['project_name'] = "FellowsYates2021"
seqdata['publication_year'] = 2021

(seqdata
    .iloc[:, [10, 11] + list(range(0, 10))]
    .to_csv(snakemake.output[0], sep="\t", index=False)
)

