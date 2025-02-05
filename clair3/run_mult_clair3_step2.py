import os
import subprocess

# Diretórios de entrada e saída
# Diretórios de entrada e saída
INPUT_DIR = "/path/to/input/directory"
OUTPUT_DIR = "/path/to/output/directory"
THREADS = "20"
MODEL_PATH = "/path/to/model/directory"
MODEL_NAME = "model_name"
DOCKER_IMAGE = "hkubal/clair3:v1.0.10"
GENOME_PATH = "/path/to/genomes/directory"
VCF_CANDIDATE_PATH = "/path/to/vcf/directory"

# Garantir que o diretório de saída exista
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Listar todos os arquivos BAM no diretório de entrada
bam_files = [filename for filename in os.listdir(INPUT_DIR) if filename.endswith(".bam")]

for bam_file in bam_files:
    # Caminhos completos
    bam_path = os.path.join(INPUT_DIR, bam_file)
    ref_path = os.path.join(GENOME_PATH, "1KG_ONT_VIENNA_hg38.fa")  # Caminho do arquivo de referência
    
    # Criar subpasta no diretório de saída com base no nome do arquivo BAM
    output_subfolder = os.path.join(OUTPUT_DIR, os.path.splitext(bam_file)[0])
    os.makedirs(output_subfolder, exist_ok=True)
    
    # Comando Docker para rodar o Clair3
    docker_command = [
        "docker", "run", "-it",
        "-v", f"{INPUT_DIR}:{INPUT_DIR}",
        "-v", f"{LONGPHASE_PATH}:{LONGPHASE_PATH}",
        "-v", f"{OUTPUT_DIR}:{OUTPUT_DIR}",
        "-v", f"{MODEL_PATH}:{MODEL_PATH}",
        "-v", f"{GENOME_PATH}:{GENOME_PATH}",
        DOCKER_IMAGE,
        "/opt/bin/run_clair3.sh",
        f"--bam_fn={bam_path}",
        f"--ref_fn={ref_path}",
        f"--threads={THREADS}",
        "--platform=ont",
        f"--model_path={MODEL_PATH}/{MODEL_NAME}",
        f"--output={output_subfolder}",
        "--include_all_ctgs",
        "--gvcf",
        "--include_all_ctgs",
        f"--vcf_fn={VCF_CANDIDATE_PATH}"
    ]
    
    try:
        # Executa o comando Docker
        subprocess.run(docker_command, check=True)
        print(f"Processed {bam_file}")
    except subprocess.CalledProcessError as e:
        print(f"Error processing {bam_file}: {e}")

print("All BAM files processed.")

