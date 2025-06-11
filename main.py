import argparse
from pathlib import Path
from reference_preparer import ReferencePreparer

def create_directory_structure():
    """Создает структуру каталогов проекта"""
    dirs = [
        'data',
        'data/genome',
        'data/annotations',
        'data/output'
    ]
    
    for directory in dirs:
        Path(directory).mkdir(parents=True, exist_ok=True)
        print(f"[✓] Создана директория: {directory}")

def main():
    """
    Основной скрипт для инициализации проекта по дизайну зондов BRCA1/2.
    Создает структуру каталогов и запускает подготовку референсных данных.
    Этапы работы:
    1. Загрузка и подготовка референсных данных (геном и аннотации)
    2. Извлечение координат экзонов BRCA1/2
    """
    parser = argparse.ArgumentParser(description='Инициализация проекта по дизайну зондов BRCA1/2')
    parser.add_argument('--force-download', action='store_true', help='Принудительная перезагрузка референсных данных')
    parser.add_argument('--force-prep', action='store_true', help='Принудительная перегенерация файлов')
    
    args = parser.parse_args()
    
    print("[1] Создание структуры каталогов...")
    create_directory_structure()
    
    print("[2] Подготовка референсных данных...")
    preparer = ReferencePreparer(
        genome_url="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz",
        gtf_url="ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
        data_dir="data"
    )
    preparer.prepare_all(force_download=args.force_download, force_preparing=args.force_prep)
    
    print("[✓] Инициализация проекта завершена")

if __name__ == "__main__":
    main()
   
   


