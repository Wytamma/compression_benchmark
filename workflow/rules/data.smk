rule get_data:
    output:
        fq=RESULTS / "data/{group}/{name}",
        size=RESULTS / "data/{group}/{name}.size"
    log:
        LOGS / "get_data/{group}/{name}.log"
    container:
        CONTAINERS["fastq_dl"]
    conda:
        ENVS / "wget.yaml"
    resources:
        mem_mb=int(2 * GB)
    params:
        opts="-v",
        outdir=lambda wildcards, output: Path(output.fq).parent,
        path=lambda wildcards, output: samples[wildcards.name]["path"]
    shadow: "shallow"
    shell:
        """
        # check if path is a url
        if [[ {params.path} =~ ^http.* ]]; then
            wget -O {output.fq} {params.path} 2>> {log}
        # check if path exists
        elif test -f {params.path}; then
            cp {params.path} {output.fq}
        else
            fastq-dl {params.opts} -a {wildcards.name} -o {params.outdir} 2> {log}
            if [ {wildcards.group} = "nanopore" ]; then 
                gzip -dc {params.outdir}/{wildcards.name}*.fastq.gz > {output.fq} 2>> {log}
            else
                READ1={params.outdir}/{wildcards.name}_1.fastq.gz
                READ2={params.outdir}/{wildcards.name}_2.fastq.gz
                (paste <(zcat $READ1 | paste - - - - ) <(zcat $READ2 | paste - - - - ) | tr '\t' '\n') > {output.fq} 2>> {log}
            fi
        fi
        (wc -c {output.fq} | awk '{{print $1}}') > {output.size} 2>> {log}
        """
