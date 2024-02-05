
rule mash_classification:
    """
    Classifies assemblies/reads with Mash
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        TSV = ROOT / 'mash' / '{class_type}' / 'dist.tsv'
    threads: 24
    params:
        db = config['db']['mash']
    benchmark:
        BENCH / 'mash_classification.{class_type}.txt'
    shell:
        """
        ml mash/2.2;
        mash dist -d 0.9 -i -p {threads} {params.db} {input} > {output.TSV};
        """
rule mash_clean:
    """
    Remove assemblies/reads with 0 matching hashes.
    Extract taxid from refseq column.
    """
    input:
        TSV = rules.mash_classification.output.TSV
    output:
        TSV = ROOT / 'mash' / '{class_type}' / 'dist_cleaned.tsv'
    log:
        LOGS / 'mash_{class_type}.log'
    run:
        import pandas as pd
        mash_output = pd.read_csv(input.TSV, sep='\t', names=['ref_id', 'query', 'dist', 'p_val', 'hashes'])
        mash_output['matching_hashes'] = mash_output['hashes'].apply(lambda x: int(x.split('/')[0]))
        # logging.info(f"{sum(mash_output['matching_hashes'] == 0)} rows dropped because zero hashes")
        ## Remove rows with zero matching hashes
        # mash_output = mash_output[mash_output['matching_hashes'] != 0]
        mash_output['tax_id'] = mash_output['ref_id'].apply(lambda x: int(x.split('|')[2]))
        with open(output.TSV,'w') as handle:
            mash_output.to_csv(handle,sep='\t',index=False, header=False)

rule mash_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        TSV = rules.mash_clean.output.TSV
    output:
        TSV = ROOT / 'mash' / '{class_type}' / 'dist_cleaned_tax.tsv'
    params:
        db=config['db']['taxonomy']
    threads: 16
    shell:
        """
        ml load taxonkit;
        paste <(cat {input.TSV}) <(cut -f 7 {input.TSV} \
        | taxonkit lineage -j {threads} --data-dir {params.db} \
        | taxonkit reformat -j {threads} --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
        | cut -f 3,4 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
        ) > {output.TSV};
        """

rule mash_filter:
    """
    Per assembly/read AND genus/species, sum number of matching hashes.
    Per assembly/read, keep genus/species with the highest matching hashes.
    If there are multiple genera/species with equal sums of matching hashes, mark assembly/read as 'ambigious'.
    """
    input:
        TSV = rules.mash_taxonomy.output.TSV
    output:
        TSV = ROOT / 'mash' / '{class_type}' / '{class_level}' /'dist_cleaned_tax_filtered.tsv'
    params:
        level = lambda wildcards: wildcards.class_level
    run:
        import pandas as pd
        mash_input = pd.read_table(input.TSV
            ,names=['ref_id', 'query', 'dist', 'p_val', 'hashes', 'matching_hashes'
                , 'tax_id', 'genus', 'species'])
        ## Per query AND genus/species, take the sum of hashes
        mash_inputg = mash_input.groupby(['query', params.level],as_index=False).agg(sum_hash=('matching_hashes', 'sum'))
        ## Select for a query the genus with the largest sum of hashes
        mash_input_max = mash_inputg.loc[mash_inputg.groupby(['query'])['sum_hash'].apply(lambda x: x == x.max())].copy()

        ## Get index of queries that have more than one equal classification
        ambig_index = mash_input_max[mash_input_max.groupby('query')[params.level].transform('size') > 1].index
        ## Update genus of those queries to 'ambigious'
        mash_input_max.loc[ambig_index, params.level] = 'ambigious'
        ## Drop the duplicate 'ambigious'
        mash_input_max.drop_duplicates(subset='query',inplace=True)


        with open(output.TSV,'w') as handle:
            mash_input_max.to_csv(handle,sep='\t',index=False, header=False)

rule mash_output:
    input:
        canu_output = lambda wildcards: FILENAME[wildcards.class_type],
        TSV_mash = rules.mash_filter.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'

    output:
        TSV_mash = ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_mash.tsv'
    params:
        level = lambda wildcards: wildcards.class_level
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        mash_input = pd.read_table(input.TSV_mash,
                                   usecols=[0, 1],
                                   names=['query', params.level])

        organism_count = mash_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        organism_count.to_csv(output.TSV_mash,sep='\t',index=True,header=False)