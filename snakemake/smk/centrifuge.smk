CENTRIFUGE = CLASSIFIERS_LOCATIONS['centrifuge']['path']
CENTRIFUGE_DB = CLASSIFIERS_LOCATIONS['centrifuge']['db_path']


rule centrifuge_classification:
    """
    Classifies reads/assemblies using centrifuge
    """
    input:
        lambda wildcards: INPUT_READS
    output:
        TSV = ROOT / 'centrifuge' / 'classification.tsv',
        TSV_summary = ROOT / 'centrifuge' / 'classification_summary.tsv'
    threads: 64
    params:
        db = CENTRIFUGE_DB,
    shell:
        """
        {CENTRIFUGE} \
        -p {threads} \
        -k 1 \
        -q \
        -x {params.db} \
        --report-file {output.TSV_summary} \
        {input} \
        > {output.TSV}
        """


rule centrifuge_clean:
    """
    Remove unclassified hits (which will be added again in the last rule as 'no_hit'
    """
    input:
        rules.centrifuge_classification.output.TSV
    output:
        TSV = ROOT / 'centrifuge' / 'classification_cleaned.tsv'
    run:
        import pandas as pd

        centrifuge_ouput = pd.read_table(input[0],
            usecols=[0, 1, 2]
        )
        ## Remove unclassified maps
        centrifuge_ouput_clean = centrifuge_ouput[centrifuge_ouput.loc[:, 'seqID'] != 'unclassified']

        with open(output.TSV,'w') as handle:
            centrifuge_ouput_clean.to_csv(handle,sep='\t',index=False)

rule centrifuge_taxonomy:
    """
   Extract genus/species names based on taxid using Taxonkit
   """
    input:
        rules.centrifuge_clean.output.TSV
    output:
        TSV = ROOT / 'centrifuge' / 'classification_cleaned_tax.tsv'
    params:
        db = TAXONOMY_DB
    shell:
        """
        ml taxonkit/0.15.0;
        FILE=$(sed '1d' {input})
        paste <(cut -f 1,3 <<<"$FILE") <(cut -f 3 <<<"$FILE" \
        | taxonkit lineage --data-dir {params.db} \
        | taxonkit reformat --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
        | cut -f 3,4 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
        ) > {output.TSV};
        """

rule centrifuge_filter:
    """
    Adds reads that were not classified.
    Split file to either taxid with genus name and taxid with species name.
    """
    input:
        TSV_centrifuge = rules.centrifuge_taxonomy.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'

    output:
        TSV = ROOT / 'output' / '{class_level}' / 'output_centrifuge.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names=lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_centrifuge.log'
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        centrifuge_input = pd.read_table(input.TSV_centrifuge,
                                         names=['query', 'tax_id', 'genus', 'species']
        )

        organism_count = centrifuge_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            organism_count=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV,sep='\t',index=True,header=False)

