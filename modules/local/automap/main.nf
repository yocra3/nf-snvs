process AUTOMAP {
    tag "$meta.id"
    label 'process_single'
    errorStrategy 'retry'

    // Conda is not supported
    container "/mnt/tblab/yolanda/automap/automap.sif"
    
    input:
    tuple val(meta), path(vcf), path(tbi)
    val automap_assembly
    path projectDir

    output:
    tuple val(meta), path("*HomRegions*"), emit: roh_automap
    path "versions.yml"           , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    prefix = task.ext.prefix ?: "${meta.id}"
    
    if ( task.attempt == 1 )
        """
        if [[ \$(bcftools query -l ${vcf} | wc -l) -gt 1 ]]; then
		    for sample in \$(bcftools query -l ${vcf}); do

			    bcftools view -s \${sample} -O v -o \${sample}.indv.vcf ${vcf}
					
			    bash /app/AutoMap_v1.3.sh \\
			    --vcf \${sample}.indv.vcf \\
			    --out . \\
			    --genome ${automap_assembly}

			    mv \${sample}/*HomRegions* .

		    done
			
	    else
					
		    bcftools view -O v -o ${vcf}.vcf ${vcf}

		    bash /app/AutoMap_v1.3.sh \\
		    --vcf ${vcf}.vcf \\
		    --out . \\
		    --genome ${automap_assembly}

		    mv \$(bcftools query -l ${vcf})/*HomRegions* .

	    fi
        
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
	    """
    
    else
        """
        echo "Less than 10k variants (with quality)" > ${prefix}_no_automap_HomRegions.txt
		
        cat <<-END_VERSIONS > versions.yml
        "${task.process}":
            bcftools: \$(bcftools --version 2>&1 | head -n1 | sed 's/^.*bcftools //; s/ .*\$//')
        END_VERSIONS
        """
}
