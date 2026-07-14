## Overview of the folder `05-results`

This folder contains all results generated in this project, including processed data and summary tables produced by the analysis workflows.

Files that are used as input for generating figures:
- `AC_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; ancient calculus MAGs
- `BM_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; buccal mucosa MAGs
- `HP_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; hard palate MAGs
- `KG_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; keratinized gingiva MAGs
- `MC_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; modern calculus MAGs
- `PQ_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; plaque MAGs
- `RS_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; Reference HOMD genomes
- `SL_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; saliva MAGs
- `TN_cas_k_o_p.tsv` - CCTyper Cas gene hits known, orphan, putative; tongue MAGs
- `card_rgi_abx_resistance_gene_table.tsv.gz` - CARD RGI antibiotic resistance gene table
- `class_order_table.tsv` - Classes listed in order by Phylum for plotting
- `dRep_rgig7111_anaerovoracaceae_Cdb.csv` - dRep cluster table for Anaerovoracaceae genomes
- `dRep_rgig7111_anaerovoracaceae_Mdb.csv.gz` - dRep MASH table for Anaerovoracaceae genomes
- `dRep_rgig7111_anaerovoracaceae_genomeInfo.csv` - dRep genome info table for Anaerovoracaceae genomes
- `doralis_scc.raxml.support` - Desulfobulbus oralis RAxML tree with bootstrap support
- `drug_classes.tsv` - table organizing CARD drug classes
- `duul01_mad_root.nwk` - DUUL01 RAxML tree with bootstrap support rooted using the MAD technique
- `genomad_phage_imgvr_taxonomy.tsv` - Taxonomy of phages identified by geNomad, classified with IMGVR4
- `human_pydamage_table.tsv.gz` - PyDamage estimations for contigs from human calculus assemblies
- `kraken2_ci15_all_oral_samples.homd_mags_11_2025_species.tsv` - species count table for all oral samples classified by Kraken2 with HM database
- `kraken2_classification_counts_gtdb207.tsv` - Number of species per oral samples classified by Kraken2 with the GTDB r207 database
- `kraken2_classification_counts_homd.tsv` - Number of species per oral samples classified by Kraken2 with the HOMD database
- `kraken2_classification_counts_homd_hpd_csc.tsv` - Number of species per oral samples classified by Kraken2 with the HOMD database, 4 additional samples
- `kraken2_classification_counts_homd_mags.tsv` - Number of species per oral samples classified by Kraken2 with the HM database
- `kraken2_classification_counts_homd_mags_csc_hpd.tsv` - Number of species per oral samples classified by Kraken2 with the HM database, 4 additional samples
- `mag_sample_to_eager_label.tsv` - Sample label conversion table
- `methanocatella_c4_mad_root.nwk` - Phylophlan tree for Methanobrevibacter cluster 1190 rooted using the MAD technique
- `neanderthal_pydamage_table.tsv` - PyDamage estimations for contigs from Neanderthal calculus assemblies
- `phables_phage_imgvr_taxonomy.tsv` - Taxonomy of phages identified by Phables, classified with IMGVR4 
- `phage_blastn_ani_clusters.tsv` - Clustering results of all complete, high, and medium quality phages identified by geNomad and Phables
- `phylophlan_aln_methanomethylophilus_JTT_w_CSC_20251105.raxml.support` - Phylophlan tree of Methanomethylophilus with bootstrap support
- `rgig7111_09_2025_mad_root.nwk` - RAxML tree of Ca. Cryptocalculus aphthonus rooted usign the MAD technique
- `rgig7111_cctyper_table.tsv` - CCTyper Cas gene hits for Ca. Cryptocalculus aphthonus
- `rgig7111_dlt_operon_table.tsv` - dlt gene table for Ca. Cryptocalculus aphthonus
- `rgig7111_van_gene_table.tsv` - van operon gene table for Ca. Cryptocalculus aphthonus
- `vanish_blossum_table.tsv` - Differential abundance table for VANISH/BloSSUM analysis
- `vb_species_plot_order.tsv` - VANISH/BloSSUM species ordered for plotting

### `Archaea`: results and tables for *Methanocatella* and *Methanomethylophilus* analyses 

*Methanocatella*:

- `methanocatella_metadata.tsv`: *Methanocatella* MAGs metadata table
- `methanocatella_contig_list.tsv`: List of contigs assigned to the genus *Methanocatella* used to extract Pydamage results
- `methanocatella_full_tree_panaroo_gene_presence_absence.Rtab`: *Methanocatella* Panaroo gene presence/absence table
- `phylophlan_aln_methanocatella_full_tree.raxml.support`: *Methanocatella* PhyloPhlAn tree with bootstraps
- `methanocatella_cctyper_cas_table.tsv`: Cas subtypes assigned to *Methanocatella* MAGs
- `methanocatella_oralis_cctyper_cas_table.tsv`: Cas subtypes assigned to *Methanocatella* *oralis* reference genome (GCF_001639275.1)  
- `methanocatella_dbcan3_cazyme_table.tsv`: Carbohydrate active enzymes (CAZymes) identified in *Methanocatella* MAGs
- `cazyme_functional_categories.tsv`: Table assigning specific CAZymes identified in archaea MAGs to broader functional categories
- `dbcan3_cazyme_table_methanocatella_oralis.tsv`: Prophages identified in *Methanocatella* MAGs
- `methanocatella_card_table.tsv`: Antibiotic resistance genes identified in *Methanocatella* MAGs
- `methanocatella_prophages.tsv`: Prophages identified in *Methanocatella* MAGs
- `methanocatella_methanogenesis_genes_presence_absence.tsv`: The presence/absence of genes involved in methane production across *Methanocatella* MAGs
- `SRR13263123.mapped_methanocatella_c90_q25l30_rmdup.bam`: BAM file of the El Sidrón Neanderthal (SRS7890498) reads mapped to a Ca. *Methanocatella* *cohabitans* cluster representative (ERS6205388.bin.r.4)
- `SRR13263123.mapped_m_oralis_q25l30_rmdup.bam`: BAM file of the El Sidrón Neanderthal (SRS7890498) reads mapped to a *Methanocatella* *oralis* genome (GCF_001639275.1)
- `cas_IIIA_percent_identity_matrix.csv`: Percent identity matrix for ten references containing Cas subtype IIIA operons
- `methanocatella_casIIIA_MAFFT_matrix.csv`: Percent identity matrix of *Methanocatella* Cas IIIA MAFFT (nucleotide) alignment 
- `methanocatella_cas_IIIA_nodes_structures.csv`: Table grouping *Methanocatella* Cas IIIA operons into structural subtypes 
- `methanocatella_casIIIA_MUSCLE_matrix.csv`: Percent identity matrix of *Methanocatella* Cas IIIA MUSCLE (protein) alignment
- `methanocatella_dRep_ANI_Mdb.csv`: Similarity between all *Methanocatella* MAGs (dRep output)
- `methanocatella_damage_curves.tsv`: PyDamage estimations for contigs from *Methanocatella* assemblies
- `methanogenesis_genes.csv`: KEGG genes used to determine methanogenesis operons
- `methanogenesis_locus_tags_and_panaroo_group.tsv`: correlating named methanogenesis genes with panaroo gene IDs

*Methanocatella* *massiliense* (pr_cl_1187):
- `methanocatella_c1187_scc.raxml.support`: *M. massiliense* (pr_cl_1187) non-recombinant single copy core gene tree with bootstraps 
- `methanocatella_c1187_mad_root.nwk`: MAD-rooted *M. massiliense* (pr_cl_1187) non-recombinant single copy core gene tree newick file
- `fastBAPS_methanocatella_cluster1187_snp_counts_and_clusters.tsv`: *M. massiliense* (pr_cl_1187) MAGs fastBAPS cluster identifications and SNP counts
- `methanocatella_c1187_non_recomb_scc_gene_list.tsv`: List of *M. massiliense* (pr_cl_1187) non-recombinant single copy core genes 
- `methanocatella_c1187_panaroo_gene_presence_absence.Rtab`: *M. massiliense* (pr_cl_1187) Panaroo gene presence/absence table

*Methanocatella* *oralis* (pr_cl_1188):
- `methanocatella_c1188_scc.raxml.support`: *M. oralis* (pr_cl_1188) non-recombinant single copy core gene tree with bootstraps 
- `methanocatella_c1188_mad_root.nwk`: MAD-rooted *M. oralis* (pr_cl_1188) non-recombinant single copy core gene tree newick file
- `fastBAPS_methanocatella_cluster1188_snp_counts_and_clusters.tsv`: *M. oralis* (pr_cl_1188) MAGs fastBAPS cluster identifications and SNP counts
- `methanocatella_C1188_non_recomb_scc_gene_list.tsv`: List of *M. oralis* (pr_cl_1188) non-recombinant single copy core genes 
- `methanocatella_c1188_panaroo_gene_presence_absence.Rtab`: *M. oralis* (pr_cl_1188) Panaroo gene presence/absence table

Ca. *Methanocatella* *senecta* (pr_cl_1189):
- `methanocatella_c1189_scc.raxml.support`: Ca. *M. senecta* (pr_cl_1189) non-recombinant single copy core gene tree with bootstraps 
- `methanocatella_c1189_mad_root.nwk`: MAD-rooted Ca. *M. senecta* (pr_cl_1189) non-recombinant single copy core gene tree newick file
- `fastBAPS_methanocatella_cluster1189_snp_counts_and_clusters.tsv`: Ca. *M. senecta* (pr_cl_1189) MAGs fastBAPS cluster identifications and SNP counts
- `methanocatella_c1189_non_recomb_scc_gene_list.tsv`: List of Ca. *M. senecta* (pr_cl_1189) non-recombinant single copy core genes 
- `methanocatella_c1189_panaroo_gene_presence_absence.Rtab`: Ca. *M. senecta* (pr_cl_1189) Panaroo gene presence/absence table

Ca. *Methanocatella* *cohabitans* (pr_cl_1190):
- `methanocatella_c1190_scc.raxml.support`: Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core gene tree with bootstraps 
- `methanocatella_c1190_mad_root.nwk`: MAD-rooted Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core gene tree newick file
- `fastBAPS_methanocatella_cluster1190_snp_counts_and_clusters.tsv`: Ca. *M. cohabitans* (pr_cl_1190) MAGs fastBAPS cluster identifications and SNP counts
- `methanocatella_c1190_non_recomb_scc_gene_list.tsv`: List of Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core genes 
- `methanocatella_c1190_panaroo_gene_presence_absence.Rtab`: Ca. *M. cohabitans* (pr_cl_1190) Panaroo gene presence/absence table
- `methanocatella_c1190_no_outgroup.caster.nwk`: Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core gene tree CASTER-site newick file with no outgroup
- `methanocatella_c1190_w_outgroup.caster.nwk`: Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core gene tree CASTER-site newick file with *M. oralis* as the outgroup
- `methanocatella_c1190_w_outgroup_scc.raxml.support`: Ca. *M. cohabitans* (pr_cl_1190) non-recombinant single copy core gene tree RAxML file with *M. oralis* as the outgroup

*Methanobacteriaceae*:
- `methanobacteriaceae_tree_metadata.tsv`: NCBI genomes and *Methanocatella* MAGs used in *Methanobacteriaceae* tree 
- `phylophlan_aln_methanobacteriaceae_SI_41.raxml.support`: *Methanobacteriaceae* PhyloPhlAn tree with bootstraps

*Methanomethylophilus*:

- `methanomethylophilus_metadata.tsv`: *Methanomethylophilus* MAGs metadata table
- `methanomethylophilus_contig_list.tsv`: List of contigs assigned to the genus *Methanomethylophilus* used to extract Pydamage results
- `methanomethylophilus_full_tree_panaroo_gene_presence_absence.Rtab`: *Methanomethylophilus* Panaroo gene presence/absence table
- `phylophlan_aln_methanomethylophilus_full_tree.raxml.support`: *Methanomethylophilus* PhyloPhlAn tree with bootstraps
- `methanomethylophilus_cctyper_cas_table.tsv`: Cas subtypes assigned to *Methanomethylophilus* MAGs
- `methanomethylophilus_alvi_cctyper_cas.tsv`: Cas subtypes assigned to *Methanomethylophilus alvi* reference genome (GCF_003711245.1)
- `methanomethylophilus_dbcan3_cazyme_table.tsv`: Carbohydrate active enzymes (CAZymes) identified in *Methanomethylophilus* MAGs
- `methanomethylophilus_prophages.tsv`: Prophages identified in *Methanomethylophilus* MAGs
- `methanomethylophilus_methanogenesis_genes_presence_absence.tsv`: The presence/absence of genes involved in methane production across *Methanomethylophilus* MAGs
- `methanomethylophilus_casIIA_MAFFT_matrix.csv`: Percent identity matrix of *Methanomethylophilus* Cas IIA MAFFT (nucleotide) alignment 
- `methanomethylophilus_cas_IIA_nodes_structures.csv`: Table grouping *Methanomethylophilus* Cas IIA operons into structural subtypes 
- `methanomethylophilus_casIIA_MUSCLE_matrix.csv`: Percent identity matrix of *Methanomethylophilus* Cas IIA MUSCLE (protein) alignment
- `methanomethylophilus_damage_curves.tsv`: PyDamage estimations for contigs from *Methanomethylophilus* assemblies

Ca. *Methanomethylophilus fidelis* (pr_cl_3192):
- `methanomethylophilus_c3192_scc.raxml.support`: Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core gene tree with bootstraps 
- `methanomethylophilus_c3192_mad_root.nwk`: MAD-rooted Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core gene tree newick file
- `fastBAPS_methanomethylophilus_cluster3192_snp_counts_and_clusters.tsv`: Ca. *M. fidelis* (pr_cl_3192) MAGs fastBAPS cluster identifications and SNP counts
- `methanomethylophilus_c3192_non_recomb_scc_gene_list.tsv`: List of Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core genes 
- `methanomethylophilus_c3192_panaroo_gene_presence_absence.Rtab`: Ca. *M. fidelis* (pr_cl_3192) Panaroo gene presence/absence table
- `methanomethylophilus_c3192_no_outgroup.caster.nwk`: Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core gene tree CASTER-site newick file with no outgroup
- `methanomethylophilus_c3192_w_outgroup_scc.raxml.support`: Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core gene tree RAxML file with *M*. pr_cl_3194 as the outgroup
- `methanomethylophilus_c3192_w_outgroup.caster.nwk`: Ca. *M. fidelis* (pr_cl_3192) non-recombinant single copy core gene tree CASTER-site newick file with *M*. pr_cl_3194 as the outgroup

*nimI*:
- `nimI_GTR.raxml.support`: *nimI* nucleotide alignment RAxML with bootstraps
- `nimI_muscle.raxml.support`: *nimI* protein alignment RAxML with bootstraps
- `nimI_CARD_table_w_taxonomy.tsv`: CARD *nimI* hits with associated taxonomy to annotate trees
