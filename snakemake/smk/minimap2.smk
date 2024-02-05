rule minimap2_classification:
    """
    Classifies reads/assemblies with Minimap2
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        TSV = ROOT / 'minimap2' / '{class_type}' / 'report.tsv'
    threads: 96
    params:
        db = config['db']['minimap2']
    shell:
        """
        ml minimap2/2.24;
        minimap2 -x map-ont {params.db} {input} -t {threads} > {output.TSV};
        """

rule minimap2_clean:
    """
    Only keep mappings with max quality per query.
    Extract taxid from refseq column.
    """
    input:
        TSV = rules.minimap2_classification.output.TSV
    output:
        TSV = ROOT / 'minimap2' / '{class_type}' / 'report_cleaned.tsv'
    run:
        import pandas as pd
        minimap2_output = pd.read_table(input.TSV
                                        ,names=['query', 'target', 't_length', 'match_bases', 'map_length'
                                                , 'map_quality']
                                        ,usecols=[0, 5, 6, 9, 10, 11])

        ## Only keep max values of 'map quality' per query
        minimap2_out_max = minimap2_output.loc[minimap2_output.groupby('query')['map_quality'].apply(lambda x: x == x.max())].copy()
        ## Extract taxid column
        minimap2_out_max['tax'] = minimap2_out_max['target'].apply(lambda x: int(x.split('|')[2]))

        with open(output.TSV,'w') as handle:
            minimap2_out_max.to_csv(handle,sep='\t',index=False, header=False)

rule minimap2_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        TSV = rules.minimap2_clean.output.TSV
    output:
        TSV = ROOT / 'minimap2' / '{class_type}' / 'report_cleaned_tax.tsv'
    params:
        db = config['db']['taxonomy']
    shell:
        """
        ml load taxonkit;
        paste <(cat {input.TSV}) <(cut -f 7 {input.TSV} \
        | taxonkit lineage --data-dir {params.db} \
        | taxonkit reformat --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified"\
        | cut -f 3,4 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
        ) > {output.TSV};
        """
rule minimap2_filter:
    """
    Per assembly/read, count occurences of genus/species.
    Per assembly/read, keep genus/species with the most occurences.
    If there are multiple genera/species with equal occurences, mark assembly/read as 'ambigious'.
    """
    input:
        TSV = rules.minimap2_taxonomy.output.TSV
    output:
        TSV = ROOT / 'minimap2' / '{class_type}' / '{class_level}' / 'report_cleaned_tax_filtered.tsv'
    params:
        level = lambda wildcards: wildcards.class_level
    run:
        import pandas as pd
        minimap2_input = pd.read_table(input.TSV
            ,names=['query', 'target', 't_length', 'match_bases', 'map_length'
                , 'map_quality', 'tax', 'genus', 'species']
        )

        ## Per query, count number of genera/species
        minimap2_inputg = minimap2_input.groupby(['query', params.level],as_index=False).agg(count=('map_quality', 'size')
            ,max_quality=('map_quality', 'max')
            ,sum_quality=('map_quality', 'sum')
        )
        ## Per query, keep genus/species with the highest count
        minimap2_input_max = minimap2_inputg.loc[
            minimap2_inputg.groupby(['query'])['count'].apply(lambda x: x == x.max())].copy()

        ## Get index of queries that have more than one equal classification
        ambig_index = minimap2_input_max[minimap2_input_max.groupby('query')[params.level].transform('size') > 1].index
        ## Update genus/species of those queries to 'ambigious'
        minimap2_input_max.loc[ambig_index, params.level] = 'ambigious'
        ## Drop the duplicate 'ambigious'
        minimap2_input_max.drop_duplicates(subset='query',inplace=True)

        with open(output.TSV,'w') as handle:
            minimap2_input_max.to_csv(handle,sep='\t',index=False, header=False)
    
rule minimap2_output:
    input:
        canu_output = lambda wildcards: FILENAME[wildcards.class_type],
        TSV_minimap2 = rules.minimap2_filter.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'

    output:
        TSV_minimap2 = ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_minimap2.tsv'
    params:
        level = lambda wildcards: wildcards.class_level
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        minimap2_input = pd.read_table(input.TSV_minimap2,
                                        usecols=[0, 1],
                                        names=['query', params.level])

        organism_count = minimap2_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        organism_count.to_csv(output.TSV_minimap2,sep='\t',index=True,header=False)

