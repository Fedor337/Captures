import argparse
from reference_preparer import ReferencePreparer
from probe_generator import ProbeGenerator
from temperature_filtered import main as filter_tm
from GC_filtered import main as filter_gc
from structure_filtered import main as filter_structure
from repetitions_filtered import main as filter_repeats
from pathlib import Path

def main():
    """
    Командный конвейер для дизайна гибридизационных зондов BRCA1/2.

    Этапы работы:
    1. Загрузка и подготовка референсных данных (геном и аннотации)
    2. Извлечение координат экзонов BRCA1/2
    3. Генерация перекрывающихся зондов заданной длины
    4. Применение фильтров: GC-состав, Tm, повторы, вторичная структура
    """
    
if __name__ == "__main__":
    main()
