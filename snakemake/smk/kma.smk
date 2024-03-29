KMA = TOOLS_LOCATIONS['kma']['path']
KMA_DB = TOOLS_LOCATIONS['kma']['db_path']

rule kma_mapping:
    input:
        lambda wildcards: INPUT_READS
    output:
        frag = ROOT / 'kma' / 'kma.frag.gz',
        mapstat = ROOT / 'kma' / 'kma.mapstat',
        res = ROOT / 'kma' / 'kma.res'
    params:
        DB = KMA_DB,
        prefix = lambda wildcards: ROOT / 'kma' / 'kma',
        MP = 20,
        MRS = 0.0,
        BC = 0.7,
        TMP = f"{ROOT / 'kma'}/"
    threads: 24
    shell:
        """
        mkdir -p {params.TMP};
        {KMA} \
        -i {input} \
        -o {params.prefix} \
        -t_db {params.DB} \
        -t {threads} \
        -mrs {params.MRS} \
        -bcNano \
        -bc {params.BC} \
        -ef \
        -a \
        -mem_mode \
        -tmp {params.TMP} \
        -1t1 \
        -matrix \
        -verbose \
        """

rule unpack_frag_gz:
    input:
        rules.kma_mapping.output.frag
    output:
        frag = ROOT / 'kma' / 'kma.frag'
    shell:
        """
        gzip -dk {input}
        """

rule kma_clean:
    """
    From the frag file, extract all aligned pairs (template, read).
    Counting each pair is the same as the 'readCountAln' column in the .mapstat file.
    One has the choice to either count the mapped reads or the aligned reads.
    I choice for the aligned reads (as some mapped reads are rejected because of a low aligning score).
    Additionally, CCMetagen filters based on the consensus which is constructed from aligned reads.
    """
    input:
        rules.unpack_frag_gz.output.frag
    output:
        TSV = ROOT / 'kma' / 'kma_cleaned.tsv'
    run:
        import pandas as pd

        kma = pd.DataFrame()

        frags = pd.read_table(input[0],
            header=None,
            usecols=[1, 2, 3, 4, 5, 6],
            names=['eq_well_mapping', 'al_score', 'start_al', 'end_al', 'template', 'query'])

        kma['template'] = frags.loc[:, ['template']]
        kma['query'] = frags.loc[:, 'query'].apply(lambda x: x.split()[0])


        kma['taxid'] = kma['template'].apply(lambda x: int(x.split('|')[2].split()[0]))
        kma.to_csv(output.TSV, sep='\t', index=False, header=False)



rule kma_taxonomy:
    """
    Extract genus/species based on taxid using Taxonkit
    """
    input:
        TSV = rules.kma_clean.output.TSV
    output:
        TSV = ROOT / 'kma' / 'kma_cleaned_tax.tsv'
    params:
        db = TAXONOMY_DB
    threads: 16
    shell:
        """
        paste <(cat {input.TSV}) <(cut -f 3 {input.TSV} \
        | {TAXONKIT} lineage -j {threads} --data-dir {params.db} \
        | {TAXONKIT} reformat -j {threads} --data-dir {params.db} --format "{{g}}\t{{s}}" --miss-rank-repl "unclassified" \
        | cut -f 3,4 \
        | awk 'BEGIN{{FS="\t"; OFS="\t"}} {{if ($1 == "unclassified" && $2 != "unclassified") {{print "missing_genus", $2}} else {{print}}}}' \
        ) > {output.TSV};
        """

rule kma_output:
    input:
        TSV_KMA = rules.kma_taxonomy.output.TSV,
        read_count= ROOT / 'input' / 'HQ.stats'

    output:
        TSV_KMA = ROOT / 'output' / '{class_level}' / 'output_kma.tsv'
    params:
        level=lambda wildcards: wildcards.class_level,
        duplicate_names= lambda wildcards: ROOT / 'output' / wildcards.class_level / 'corrected_duplicates_kma.log'
    run:
        import pandas as pd

        with open(input.read_count,'r') as f:
            next(f)
            read_count = int(f.read().split('\t')[3])

        kma_input = pd.read_table(input.TSV_KMA,
            names=['template', 'query', 'tax_id', 'genus', 'species'])

        organism_count = kma_input[params.level].value_counts()
        organism_count['no_hit'] = read_count - organism_count.sum()

        if any(organism_count.index.duplicated(keep=False)):
            organism_count[organism_count.index.duplicated()].to_csv(params.duplicate_names, sep='\t',index=True,header=False)
            organism_count=organism_count.groupby(organism_count.index,sort=False).sum()

        # Save output
        organism_count.to_csv(output.TSV_KMA,sep='\t', index=True, header=False)

if 'ccmetagen' in USE_CLASSIFIERS:
    include: 'ccmetagen.smk'

