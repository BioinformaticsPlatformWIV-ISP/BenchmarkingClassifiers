KAIJU = CLASSIFIERS_LOCATIONS['kaiju']['path']
KAIJU_DB = CLASSIFIERS_LOCATIONS['kaiju']['db_path']


rule kaiju_classification:
    """
    Classifies reads/assemblies using Kaiju
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        TSV = ROOT / 'kaiju' / '{class_type}' / 'classification.tsv'
    params:
        db = KAIJU_DB,
        kaiju_tax = Path(TAXONOMY_DB) / 'nodes.dmp'
    threads: 16
    shell:
        """
        {KAIJU} \
        -t {params.kaiju_tax} \
        -f {params.db} \
        -i {input} \
        -v \
        -o {output} \
        -z {threads} \
        """

rule kaiju_clean:
    """
    Remove all reads/assemblies which are marked 'unclassified'
    """
    input:
        TSV = rules.kaiju_classification.output.TSV
    output:
        TSV_classified = ROOT / 'kaiju' / '{class_type}' / 'classified.tsv',
        TSV_unclassified= ROOT / 'kaiju' / '{class_type}' / 'unclassified.tsv',

    run:
        import pandas as pd

        classified = {'name': [],
                      'taxonid': [],
                      'score': [],
                      'matches_tax': []}
        unclassified = {'name': []}
        with open(input.TSV, 'r') as f:
            for line in f:
                if line[0] == 'C':
                    columns = line.split('\t')
                    classified['name'].append(columns[1])
                    classified['taxonid'].append(columns[2])
                    classified['score'].append(columns[3])
                    classified['matches_tax'].append(columns[4][:-1])
                else:
                    columns = line.split("\t")
                    unclassified['name'].append(columns[1])

        pd.DataFrame.from_dict(classified).to_csv(output.TSV_classified,sep='\t',index=False, header=False)
        pd.DataFrame.from_dict(unclassified).to_csv(output.TSV_unclassified,sep='\t',index=False, header=False)

rule kaiju_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        TSV_classified = rules.kaiju_clean.output.TSV_classified
    output:
        TSV_classified_tax = ROOT / 'kaiju' / '{class_type}' / 'classified_tax.tsv'
    params:
        TAXONOMY = TAXONOMY_DB
    shell:
        """
        ml taxonkit/0.15.0;
        export TAXONKIT_DB={params.TAXONOMY}
        cut -f 1,2 {input.TSV_classified} | \
        taxonkit lineage -r -L -i 2 | taxonkit reformat --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" -I 2 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($4 == "unclassified" && $5 != "unclassified") {{print $1, $2, $3, "missing_genus", $5}} else {{print}}}}' \
        > {output.TSV_classified_tax};
        """


rule kaiju_output:
    input:
        TSV_unclassified = rules.kaiju_clean.output.TSV_unclassified,
        kaiju_taxonomy = rules.kaiju_taxonomy.output.TSV_classified_tax
    output:
        TSV_kaiju = ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_kaiju.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_type / wildcards.class_level / 'corrected_duplicates_kaiju.log'
    run:
        import pandas as pd

        kaiju_classified = pd.read_table(input.kaiju_taxonomy,
                           names=['query', 'tax_id', 'rank', 'genus', 'species']
                           )
        kaiju_unclassified = pd.read_table(input.TSV_unclassified,
                             names=['query']
                             )

        kaiju_unclassified[['genus', 'species']] = 'no_hit'

        kaiju_all = pd.concat([kaiju_classified, kaiju_unclassified], ignore_index=True)

        organism_count = kaiju_all[params.level].value_counts()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            organism_count=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_kaiju,sep='\t',index=True,header=False)
