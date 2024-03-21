METAPHLAN_ENV = TOOLS_LOCATIONS['metaphlan']['path_env']
METAPHLAN_DB = TOOLS_LOCATIONS['metaphlan']['db_path']
BOWTIE2 = TOOLS_LOCATIONS['bowtie2']['path']


rule metaphlan_classification:
    """
    Classifies reads/assemblies using metahplan
    """
    input:
        lambda wildcards: INPUT_READS
    output:
        TSV = ROOT / 'metaphlan' / 'classification.tsv',
        bowtie_out = ROOT / 'metaphlan' / 'classification.bow.txt'
    params:
        db = METAPHLAN_DB,
        env = str(Path(METAPHLAN_ENV) / 'bin' / 'activate'),
        alignement_type = 'sensitive-local'
    shell:
        """
        export PATH=""$(dirname "{BOWTIE2}")":$PATH"
        . {params.env}
        metaphlan \
        --input_type fastq \
        --bowtie2db {params.db} \
        --nproc 16 \
        --bt2_ps {params.alignement_type} \
        --bowtie2out {output.bowtie_out} \
        --unknown_estimation \
        {input} \
        {output.TSV};
        deactivate;
        """

rule metaphlan_clean:
    """
    Cleans output of metaphlan.
    By default metaphlan outputs four columns: clade_name, NCBI_tax_id, relative_abundance, additional_species.
    This rules transform the output to read_id, taxid, genus name and species name
    """

    input:
        rules.metaphlan_classification.output.TSV
    output:
        TSV = ROOT / 'metaphlan' / 'classification_cleaned.tsv'
    threads: 8
    run:
        import pandas as pd
        import re


        species_taxids = []
        abundances = []
        with open(input[0], 'r') as f:
            for line in f:
                if line[0] == '#':
                    continue
                if line.split('\t')[0] == "UNKNOWN":
                    # do not add unknown, they will be added in the last rule
                    continue
                if re.search(r'\|s__',line):
                    _, lineage, abundance, _ = line.split('\t')
                    lca = int(lineage.split('|')[-1])
                    species_taxids.append(lca)
                    abundances.append(f'{float(abundance):.5f}')

        pd.Series(data=abundances,index=species_taxids).to_csv(output.TSV, sep='\t', header=False)

rule metaphlan_taxonomy:
    input:
        cleaned_abundance = rules.metaphlan_clean.output.TSV
    output:
        cleaned_abundance_tax = ROOT / 'metaphlan' / 'classification_cleaned_tax.tsv'
    params:
        db=TAXONOMY_DB
    shell:
        """
        paste <(cut -f 1 {input.cleaned_abundance} \
               | {TAXONKIT} lineage --data-dir {params.db} \
               | {TAXONKIT} reformat --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
               | cut -f 3,4 \
               | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
               ) \
              <(cut -f 2 {input.cleaned_abundance}) \
        > {output.cleaned_abundance_tax};
        """

rule metaphlan_output:
    input:
        cleaned_abundance_tax = rules.metaphlan_taxonomy.output.cleaned_abundance_tax,
        classification = rules.metaphlan_classification.output.TSV

    output:
        TSV_metaphlan = ROOT / 'output' / '{class_level}' / 'output_metaphlan.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_metaphlan.log'
    run:
        import pandas as pd

        with open(input.classification, 'r') as f:
            for line in f:
                if line.split('\t')[0] == "UNKNOWN":
                    unknown_portion = float(line.strip().split('\t')[-1])

        metaphlan_input = pd.read_table(input.cleaned_abundance_tax,
                                        names=['genus', 'species', 'abundance'])

        organism_count = pd.Series(metaphlan_input['abundance'].values / 100, index=metaphlan_input[params.level])
        organism_count['no_hit'] = unknown_portion / 100

        # Sum for the same organisms. This is important with genera as there can be multiple species with the same genus.
        organism_count = organism_count.groupby(organism_count.index).sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            braken2_input=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_metaphlan, sep='\t', index=True, header=False)

