CCMETAGEN_ENV = CLASSIFIERS_LOCATIONS['ccmetagen']['path_env']


rule ccmetagan_filtering:
    input:
        res = rules.kma_mapping.output.res
    output:
        results = ROOT / 'ccmetagen' / 'results.ccm.csv',
    params:
        output_prefix = lambda wildcards: ROOT / 'ccmetagen' / 'results'
    shell:
        """
        . {CCMETAGEN_ENV}/bin/activate
        mkdir -p $(dirname {params.output_prefix})
        CCMetagen.py -i {input.res} \
                     -o {params.output_prefix} \
                     -r RefSeq \
                     -m text;
        deactivate;
        """

rule ccmetagen_clean:
    """
    Remove all unnecessary columns from ccmetagen's output.
    Extract from kma's mapstat the number of aligned reads to the templates and merge with ccmetagen's output.   
    """
    input:
        results = rules.ccmetagan_filtering.output.results,
        mapstat = rules.kma_mapping.output.mapstat
    output:
        results_cleaned = ROOT / 'ccmetagen' / 'results_cleaned.tsv',
    run:
        import pandas as pd
        results = pd.read_table(input.results, sep=',', usecols=[0, 11])
        kma_mapstat = pd.read_table(input.mapstat, skiprows=6, usecols=[0, 1, 2, 13, 14])
        merged = results.merge(kma_mapstat[
            ["# refSequence", 'readCountAln']],how='left',left_on='Closest_match',right_on='# refSequence').drop(columns='# refSequence')

        merged.to_csv(output.results_cleaned, sep='\t', index=False, header=False)

rule ccmetagen_taxonomy:
    input:
        results_cleaned = rules.ccmetagen_clean.output.results_cleaned
    output:
        results_cleaned_tax = ROOT / 'ccmetagen' / 'results_cleaned_tax.tsv',
    params:
        db = TAXONOMY_DB
    shell:
        """
        ml taxonkit/0.15.0;
        paste <(cut -f 1,2 {input}) \
              <(cut -f 2 {input} \
                | taxonkit reformat -I 1 --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
                | cut -f 2,3 \
                | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}') \
              <(cut -f 3 {input}) \
              > {output.results_cleaned_tax};
        """

rule ccmetagen_output:
    input:
        TSV_ccmetagen = rules.ccmetagen_taxonomy.output.results_cleaned_tax,
        read_count = ROOT / 'input' / 'HQ.stats'
    output:
        TSV_ccmetagen = ROOT / 'output' / '{class_level}' / 'output_ccmetagen.tsv'
    params:
        level = lambda wildcards: wildcards.class_level,
        duplicate_names=lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_ccmetagen.log'
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        ccmetagen_input = pd.read_table(input.TSV_ccmetagen,
                        names=['reference', 'tax_id', 'genus', 'species', 'count']
                    )

        organism_count = pd.Series(ccmetagen_input['count'].values, index=ccmetagen_input[params.level])
        organism_count = organism_count.groupby(level=0).sum()
        organism_count['no_hit'] = read_count - organism_count.sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            organism_count=organism_count.groupby(organism_count.index,sort=False).sum()

        organism_count.to_csv(output.TSV_ccmetagen,sep='\t',index=True,header=False)

