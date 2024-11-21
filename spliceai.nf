process SpliceAI {
    // Definir el comando para ejecutar SpliceAI aquí
    input:
    path vcf_file

    output:
    path "output_*.txt"

    script:
    """
    spliceai --input ${vcf_file} --output output_results.txt
    """
}
