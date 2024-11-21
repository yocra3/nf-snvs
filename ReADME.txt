## Requisitos

Este pipeline requiere que `SpliceAI` esté disponible en el entorno. Puedes instalarlo utilizando pip o Conda.

### Instalación con pip

Si usas pip, ejecuta:

```bash
pip install -r requirements.txt


### Realizar Pruebas

Antes de hacer el push al repositorio, prueba el flujo de trabajo localmente para asegurarte de que todo funcione como se espera. Puedes ejecutar tu flujo de trabajo con el siguiente comando:

```bash
nextflow run main.nf --vcf_file /ruta/al/archivo.vcf
