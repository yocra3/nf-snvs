include { SpliceAI } from './modules/spliceai/spliceai.nf'

// En tu flujo de trabajo principal
workflow {
    SpliceAI(vcf_file: 'ruta/al/archivo.vcf')
}

