rule mmseqs2_classification:
    """
    Classifies reads/assemblies using mmseqs2
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        TSV = ROOT / 'mmseqs2' / '{class_type}' / 'classification_lca.tsv',
    threads: 64
    params:
        db=config['db']['mmseqs2'],
        prefix = lambda wildcards: ROOT / 'mmseqs2' / wildcards.class_type / 'classification',
        tmp = lambda wildcards: ROOT / 'mmseqs2' / wildcards.class_type / 'tmp'
    shell:
        """
        ml mmseqs2;
        mmseqs easy-taxonomy \
        --tax-lineage 1 \
        {input} \
        {params.db} \
        {params.prefix} \
        {params.tmp} \
        --threads {threads} \
        """

# rule mmseqs2_clean:
#     """
#     Extract genus and species for the reads that have a hit.
#     Reads with a classification level above genus or species will be marked 'unclassified'
#     There are two types of proteins without a hit: those that get prefiltered (so not even searched with) and those
#     that returned nothing on a hit. The latter shows up in this tsv, so they need to be filtered out.
#     """
#     input:
#         TSV_mmseqs2 = rules.mmseqs2_classification.output.TSV
#     output:
#         TSV = ROOT / 'mmseqs2' / '{class_type}' / 'classification_cleaned.tsv'
#     run:
#         import pandas as pd
#         import re
#
#         df = pd.read_table(input.TSV_mmseqs2, header=None)
#
#         # filter out proteins that got searched with but returned no hit
#         df = df[~df.iloc[:, -1].isnull()]
#
#         def organism_name(row, rank):
#             """
#             Last column of df is formatted as:
#             '-_cellular organisms;d_Bacteria;-_Terrabacteria group;p_Actinobacteria;c_Actinomycetia'
#             return the name of genus or species if 's_' or 'g_' is present, respectively. Otherwise, return
#             'unclassified'
#             """
#             short_rank = {'species': 's_',
#                           'genus': 'g_'}
#             m = re.search(short_rank[rank] + '(.*?);',row)
#             if m:
#                 return m.group(1)
#             else:
#                 return 'unclassified'
#
#
#         tsv = {'read': list(df.iloc[:, 0]),
#          'taxid': list(df.iloc[:, 1]),
#          'genus': list(df.iloc[:, -1].apply(lambda x: organism_name(x,'genus'))),
#          'species': list(df.iloc[:, -1].apply(lambda x: organism_name(x,'species')))}
#
#         pd.DataFrame.from_dict(tsv).to_csv(output.TSV, header=False, index=False, sep='\t')

rule mmseqs2_tax:
    input:
        TSV_mmseqs2=rules.mmseqs2_classification.output.TSV
    output:
        TSV = ROOT / 'mmseqs2' / '{class_type}' / 'classification_tax.tsv'
    params:
        TAXONOMY=config['db']['taxonomy']
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
        canu_output = lambda wildcards: FILENAME[wildcards.class_type],
        TSV_mmseqs = rules.mmseqs2_tax.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'

    output:
        TSV = ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_mmseqs2.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_type / wildcards.class_level / 'corrected_duplicates_mmseqs2.log'
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

