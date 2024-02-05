rule desamba_classification:
    """
    Classifies reads using deSAMBA
    """
    input:
        lambda wildcards: FILENAME[wildcards.class_type]
    output:
        SAM = ROOT / 'desamba' / '{class_type}' / 'classification.sam',
    params:
        db = config['db']['desamba'],
    threads: 32
    shell:
        """
        cd /scratch/alvanuffelen/hackaton/SOFTWARE/deSAMBA;
        ./bin/deSAMBA classify -t {threads} {params.db} {input} -o {output.SAM};
        """

rule desamba_classification_analysis:
    """
    Analyses the classification file from deSAMBA
    """
    input:
        rules.desamba_classification.output.SAM
    output:
        REPORT = ROOT / 'desamba' / '{class_type}' / 'classification.report',
    params:
        db = config['db']['desamba'],
        db_tax=Path(config['db']['taxonomy']) / 'nodes.dmp'
    shell:
        """
        cd /scratch/alvanuffelen/hackaton/SOFTWARE/deSAMBA;
        ./bin/deSAMBA analysis ana_meta {input} {params.db_tax} > {output.REPORT}
        """


rule desamba_clean:
    """
    The output of deSAMBA is formatted the same as Kraken2. However, a big difference is that deSAMBA does not work
    with an LCA algorithm. Therefore reads are not assigned to levels higher than species because all the reference
    genomes are from species (or below). Number of reads belonging to a certain rank is the cumulative sum of reads
    from lower ranks. 
    1. Removes first and last line
    2. Cut and take '{taxid}  {number of reads}'
    3. Add rank to each line
    4. Only keep lines which are species
    5. Add genus and species name to each line. If a species does not have a genus rank, it gets replace by 'missing_genus'
    """
    input:
        rules.desamba_classification_analysis.output.REPORT
    output:
        TSV = ROOT / 'desamba' / '{class_type}' / 'classification_cleaned.tsv'
    shell:
        """
        ml taxonkit;
        sed '1d;$d' {input} | \
        cut -d ':' -f 2 | \
        sed 's/  /\t/' | \
        taxonkit lineage -L -r | \
        awk 'BEGIN{{FS="\t";OFS="\t"}}{{if ($3 == "species") print}}' | \
        taxonkit reformat -I 1 -f '{{g}}\t{{s}}' -r 'missing_genus' \
        > {output.TSV}
        """

rule desamba_output:
    input:
        canu_output = lambda wildcards: FILENAME[wildcards.class_type],
        TSV_desamba = rules.desamba_clean.output.TSV,
        read_count = ROOT / 'input' / 'HQ.stats',
        classification_report = rules.desamba_classification_analysis.output.REPORT
    output:
        TSV_desamba = ROOT / 'output' / '{class_type}' / '{class_level}' / 'output_desamba.tsv'
    params:
        level = lambda wildcards: wildcards.class_level
    run:
        import pandas as pd

        with open(input.read_count, 'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        with open(input.classification_report, 'r') as f:
            next(f)
            root_read_count = int(next(f).split('  ')[-1])

        desamba_input = pd.read_table(input.TSV_desamba,
                        names=['id', 'n_reads', 'rank', 'genus', 'species']
                    )

        assert root_read_count == desamba_input['n_reads'].sum()

        # All taxids are from species, so group by genus and sum the number of reads of all species belonging to the
        # same genus.
        organism_count = desamba_input.groupby(params.level)['n_reads'].agg('sum')
        organism_count['no_hit'] = read_count - organism_count.sum()

        organism_count.to_csv(output.TSV_desamba,
            sep='\t',index=True,header=False)

