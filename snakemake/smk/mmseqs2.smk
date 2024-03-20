MMSEQS2 = CLASSIFIERS_LOCATIONS['mmseqs2']['path']
MMSEQS2_DB = CLASSIFIERS_LOCATIONS['mmseqs2']['db_path']


rule mmseqs2_classification:
    """
    Classifies reads/assemblies using mmseqs2
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        TSV = ROOT / 'mmseqs2' / 'classification_lca.tsv',
    threads: 64
    params:
        db = MMSEQS2_DB,
        prefix = lambda wildcards: ROOT / 'mmseqs2' / 'classification',
        tmp = lambda wildcards: ROOT / 'mmseqs2' / 'tmp'
    shell:
        """
        {MMSEQS2} easy-taxonomy \
        --tax-lineage 1 \
        {input} \
        {params.db} \
        {params.prefix} \
        {params.tmp} \
        --threads {threads} \
        """

rule mmseqs2_tax:
    input:
        TSV_mmseqs2=rules.mmseqs2_classification.output.TSV
    output:
        TSV = ROOT / 'mmseqs2' / 'classification_tax.tsv'
    params:
        TAXONOMY = TAXONOMY_DB
    shell:
        """
        ml taxonkit/0.15.0;
        export TAXONKIT_DB={params.TAXONOMY}
        cut -f 1,2 {input.TSV_mmseqs2} | \
        taxonkit reformat --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" -I 2 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($3 == "unclassified" && $4 != "unclassified") {{print $1, $2, "missing_genus", $4}} else {{print}}}}' \
        > {output.TSV};
        """


rule mmseqs2_filter:
    """
    Adds reads that were not classified.
    Split file to either taxid with genus name and taxid with species name.
    """
    input:
        TSV_mmseqs = rules.mmseqs2_tax.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'
    output:
        TSV = ROOT / 'output' / '{class_level}' / 'output_mmseqs2.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_mmseqs2.log'
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        mmseqs2_input = pd.read_table(input.TSV_mmseqs,
            names=['query', 'tax_id', 'genus', 'species'],
            dtype={'query': 'object',
                   'tax_id': int,
                   'genus': 'object',
                   'species': 'object'}
        )

        organism_count = mmseqs2_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names,sep='\t',index=True,header=False)
            organism_count = organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV,sep='\t',index=True,header=False)

