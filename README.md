# bioalg

Пример команды: python align.py --seq ACG ATCG --method nw --match 1 --mismatch -1 --gap -2 --output result.txt


--method nw     # для Needleman-Wunsch

--method sw     # для Smith-Waterman

--fasta	Переключает в режим чтения из FASTA-файлов	

Пример: --fasta file1.fasta file2.fasta

--seq	Переключает в режим прямого ввода последовательностей	

Пример: --seq AGTCC ATGGC

--match	Балл за совпадение	

Пример: --match 1

--mismatch	Балл за несоответствие	

Пример: --mismatch -1

--gap	Штраф за гэп	

Пример: --gap -2

--output	Файл для сохранения результата (опционально)	

Пример: --output result.txt
