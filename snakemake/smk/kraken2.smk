KRAKEN2 = CLASSIFIERS_LOCATIONS['kraken2']['path']
KRAKEN2_DB = CLASSIFIERS_LOCATIONS['kraken2']['db_path']

rule kraken2_classification:
    """
    Classifies reads/assemblies using Kraken2
    """
    input:
        lambda wildcards: INPUT_READS
    output:
        TSV = ROOT / 'kraken2' / 'classification.tsv',
        TSV_report = ROOT / 'kraken2' / 'report.tsv'
    params:
        db = KRAKEN2_DB
    threads: 32
    shell:
        """
        {KRAKEN2} {input} --db {params.db} --output {output.TSV} --report {output.TSV_report} --threads {threads}
        """

rule kraken2_clean:
    """
    Remove all reads/assemblies which are marked 'unclassified'
    """
    input:
        TSV = rules.kraken2_classification.output.TSV
    output:
        TSV = ROOT / 'kraken2' / 'classification_cleaned.tsv'
    run:
        import pandas as pd
        kraken2_ouput = pd.read_table(input.TSV,
                                      header = None,
                                      usecols=[0,1,2,3]
                                     )
        ## Remove unclassified maps
        kraken2_ouput_clean = kraken2_ouput[kraken2_ouput.iloc[:,0] == 'C']

        with open(output.TSV,'w') as handle:
            kraken2_ouput_clean.to_csv(handle,sep='\t',index=False, header=False)

rule kraken2_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        rules.kraken2_clean.output.TSV
    output:
        TSV = ROOT / 'kraken2' / 'classification_cleaned_tax.tsv'
    params:
        db = TAXONOMY_DB
    shell:
        """
        ml load taxonkit;
        paste <(cut -f 2,3 {input}) <(cut -f 3 {input} \
        | taxonkit lineage --data-dir {params.db} \
        | taxonkit reformat --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
        | cut -f 3,4 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
        ) > {output.TSV};
        """

rule kraken2_output:
    input:
        TSV_kraken2 = rules.kraken2_taxonomy.output.TSV,
        read_count = ROOT / 'input' / 'HQ.stats'
    output:
        TSV_kraken2 = ROOT / 'output' / '{class_level}' / 'output_kraken2.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_kraken2.log'
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        kraken2_input = pd.read_table(input.TSV_kraken2,
                        names=['query', 'tax_id', 'genus', 'species']
                    )

        organism_count = kraken2_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            organism_count=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_kraken2,sep='\t',index=True,header=False)


if 'bracken' in USE_CLASSIFIERS:
    include: "bracken.smk"





