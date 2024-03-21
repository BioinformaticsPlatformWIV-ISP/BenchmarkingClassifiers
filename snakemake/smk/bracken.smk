BRACKEN = TOOLS_LOCATIONS['bracken']['path']
BRACKEN_DB = TOOLS_LOCATIONS['bracken']['db_path']
KRAKEN2 = TOOLS_LOCATIONS['kraken2']['path']


PYTHON_SCRIPTS = config['python_scripts']

rule convert_kraken2_classification:
    """
    Bracken directly used the report output of Kraken2.
    It interprets the report which was build with the older taxonomy.
    Once that report it interpreted, the tax ids cant be changed anymore. 
    For example: 2933872 belongs to the Ectobacillus genus.
    However, at the time of the Kraken2_DB creation, it belonged to the genus Bacillus.
    So Bracken will assign it to Bacillus and will also receive the tax id of Bacillus (because Bracken can work on genus level). 
    Therefore, the report of kraken should first be updated before passing to Bracken
    """
    input:
        TSV=rules.kraken2_classification.output.TSV,
    output:
        TSV=ROOT / 'bracken' / 'kraken2_classification_updated.tsv',
        report=ROOT / 'bracken' / 'kraken2_report_updated.tsv',
        taxonomy_db=ROOT / 'bracken' / 'mydb_taxonomy.tsv',
    params:
        TAXONOMY=TAXONOMY_DB,
        SCRIPTS=PYTHON_SCRIPTS,
    shell:
        """
        export TAXONKIT_DB={params.TAXONOMY}
        python3 {params.SCRIPTS}/make_ktaxonomy.py --nodes {params.TAXONOMY}/nodes.dmp --names {params.TAXONOMY}/names.dmp -o {output.taxonomy_db}
        paste <(cut -f1,2 {input.TSV}) <(cut -f 3 {input.TSV} | {TAXONKIT} lineage -Lnc | awk -F '\\t' '$1==0 {{print $1;next}} {{print $2}}') \
                <(cut -f 5 {input.TSV}) > {output.TSV}
        python3 {params.SCRIPTS}/make_kreport.py -i {output.TSV} -t {output.taxonomy_db} -o {output.report}
        """

rule bracken_classification:
    """
    Bayesian Reestimation of Abundance with KrakEN.
    """
    input:
        rules.convert_kraken2_classification.output.report
    output:
        TSV=ROOT / 'bracken' / '{class_level}' / 'classification.tsv',
        report=ROOT / 'bracken' / '{class_level}' / 'classification.report',
        stdout=ROOT / 'bracken' / '{class_level}' / 'bracken.output'
    params:
        db=BRACKEN_DB,
        read_len=1000,
        level=lambda wildcards: {'genus': 'G', 'species': 'S'}[wildcards.class_level]
    shell:
        """
        export PATH=""$(dirname "{KRAKEN2}")":$PATH"
        {BRACKEN} \
        -d {params.db} \
        -i {input} \
        -o {output.TSV} \
        -w {output.report} \
        -r {params.read_len} \
        -l {params.level} \
        > {output.stdout}
        """

rule bracken_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        rules.bracken_classification.output.TSV
    output:
        TSV=ROOT / 'bracken' / '{class_level}' / 'classification_tax.tsv'
    params:
        db=TAXONOMY_DB,
        level=lambda wildcards: 'g' if wildcards.class_level == "genus" else 's'
    shell:
        """
        # save header to new file
        head -n1 {input} > {output.TSV}
        paste <(cut -f 2 {input} | tail -n+2 \
                    | {TAXONKIT} lineage --data-dir {params.db} \
                    | {TAXONKIT} reformat --data-dir {params.db} --format "{{{params.level}}}" --miss-rank-repl "unclassified" \
                    | cut -f 3 \
                ) \
              <(cut -f 2- {input} | tail -n+2) \
              >> {output.TSV};
        """

rule bracken_output:
    input:
        TSV_bracken = rules.bracken_taxonomy.output.TSV,
        stdout = rules.bracken_classification.output.stdout,
        read_count = ROOT / 'input' / 'HQ.stats'
    params:
        duplicate_names=lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_bracken.log'
    output:
        TSV_bracken=ROOT / 'output' / '{class_level}' / 'output_bracken.tsv'
    run:
        import pandas as pd
        import re

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        # Number of reads not distributed
        with open(input.stdout, 'r') as stdout:
            full_text = stdout.read()
        distr_match = re.search(r'Reads distributed: (\d+)', full_text)
        not_distr_match = re.search(r'Reads not distributed.*: (\d+)', full_text)
        reads_distributed_theoretical = int(distr_match.group(1))
        reads_not_distributed = int(not_distr_match.group(1))

        braken2_input = pd.read_table(input.TSV_bracken,
            usecols=[0, 5],
            index_col=0,
            header=0,
            names=['none', 'count']).squeeze('columns')

        # Due to rounding errors, the number of added reads to species is wrong. These reads are therefore 'lost' and
        # should be compensated for with a new category
        # https://github.com/jenniferlu717/Bracken/issues/77
        reads_distributed_practice = pd.read_table(input.TSV_bracken,
            usecols=['added_reads']).squeeze('columns').sum()
        braken2_input['rounding_error'] = reads_distributed_theoretical - reads_distributed_practice

        braken2_input['not_distributed'] = reads_not_distributed
        braken2_input['no_hit'] = read_count - braken2_input.sum()

        if any(braken2_input.index.duplicated(keep=False)):
            braken2_input[braken2_input.index.duplicated()].to_csv(params.duplicate_names,sep='\t',index=True,header=False)
            braken2_input = braken2_input.groupby(braken2_input.index,sort=False).sum()

        braken2_input.to_csv(output.TSV_bracken,sep='\t',index=True,header=False)

