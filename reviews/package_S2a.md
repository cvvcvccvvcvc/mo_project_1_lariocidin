# S2a: калибровка семейства и проверка HGT-кандидата

Дата: 10.09.2026. Статус: S2a выполнен; расширенный поиск не начат. Это пакет для review A+.

## Вопрос этапа

Достаточен ли опубликованный набор для воспроизводимого поиска lariocidin/triculamin-like циклаз и устойчив ли исходный смешанный кластер *Photorhabdus asymbiotica* + три Burkholderiales?

## Результат

| Проверка | Наблюдение | Вывод |
|---|---|---|
| Литературный reference | 393 BGC: 139 canonical, 174 noncanonical, 80 других lasso-систем; 307 уникальных явно размеченных C и 215 уникальных прекурсоров | Canonical и noncanonical образуют две разные поисковые ветви; вхождение и уникальная последовательность не смешиваются |
| Экспериментальный статус | Только 4 якоря имеют экспериментальную поддержку: lariocidin, gelatinamin, palmamin, triculamin; остальные — опубликованные предсказания | Нельзя называть 393 записи подтверждёнными продуцентами |
| Связь с исходными 29 | 25/29 LrcC имеют точное совпадение в архиве Merrild; четыре неразмеченные C-кандидата сохранены как inference, а не авторская аннотация | Исходные 29 почти полностью вложены в более широкий опубликованный набор и не измеряют новизну сами по себе |
| Кластерно-разделённая HMM-проверка | При 70% identity/80% reciprocal coverage: 19 canonical- и 31 noncanonical-кластер; пять folds. Минимум held-cluster positive score 278.4, максимум консервативно оценённого comparator score 95.3 | Внутри опубликованного набора профили хорошо разделяют 227 positives и 80 comparators; это внутренняя калибровка, не оценка природных precision/recall |
| Исходный смешанный split C | 95.7/92 на полном дереве, 94.3/86 без коротких Actinomycetota; отсутствует во всех трёх локальных Pseudomonadota-анализах, включая `-bnni` и occupancy-mask | Топологический HGT-сигнал неустойчив и сам по себе слаб |
| Независимый от топологии признак локуса | У *P. asymbiotica* ближайшие A, B1, B2, C и D среди 14 BGC принадлежат Burkholderiales. Для B2/C/D: 64.0/62.5/71.6% identity при 99–100% query coverage; A покрывает лишь 61% | Сохраняется совместимый с переносом целого локуса сигнал, но пять сцепленных генов не являются независимыми повторностями |
| Полные NCBI-контексты | Все четыре C точно совпали с reference. У *Burkholderia* IMCC1007 IS30-transposase расположен в 312 bp от края A–D-локуса, но сам локус начинается в 739 bp от края контига; в двух других контекстах мобильных keyword-аннотаций рядом нет | Случай остаётся `provisional HGT candidate`; краевой контиг и отсутствие близких родственников не позволяют утверждать HGT |

Источники reference и методической интерпретации: [Merrild et al. 2025](https://pmc.ncbi.nlm.nih.gov/articles/PMC12070363/), [OSF CZPFE](https://osf.io/czpfe/), [IQ-TREE tutorial](https://iqtree.github.io/doc/Tutorial). Полные записи: [Burkholderia](https://www.ncbi.nlm.nih.gov/nuccore/NZ_JAPVQY010000029.1), [Paraburkholderia](https://www.ncbi.nlm.nih.gov/nuccore/NZ_AXAO01000004.1), [Photorhabdus](https://www.ncbi.nlm.nih.gov/nuccore/RBLJ01000001.1), [Ralstonia](https://www.ncbi.nlm.nih.gov/nuccore/CP006668.1).

## Принятые решения

1. Искать двумя профилями — canonical и noncanonical. Результат HMM является только sequence candidate; статус семейства требует предшественник и согласованный BGC-контекст.
2. Значения 95.31 и 278.4 использовать в S3 только как предварительные границы grey/high-priority, а не как универсальные биологические пороги: они получены на ограниченном sequence-informed наборе и ещё не проверены на полном фоне MAG-протеомов.
3. Не использовать исходный четырёхтаксонный split как самостоятельное доказательство или основной результат. Сохранить локус *P. asymbiotica* для повторной проверки после добавления близких носителей/неносителей и независимого дерева хозяев.
4. HGTphyloDetect/WAAFLE пока не запускать: сначала нужны расширенный каталог, близкие родственники и пригодные контиги.

## Предлагаемый контракт S3

- Primary denominator: 399 бактериальных MAG Maize Rhizosphere v1.0; 22 изолята и 10 архейных MAG учитывать отдельно.
- Первый вход: MGnify `protein_catalogue-100.tar.gz` (около 297 MB сжатым). 100%-каталог сохраняет все уникальные белки, а membership-таблица — вхождения; 50%-кластеризация слишком груба для окончательной классификации.
- Поиск обоими HMM; сохранять все occurrence mappings. Предварительные уровни: `<95.31` — background по текущей калибровке, `95.31–278.39` — grey, `≥278.4` — high sequence support. Границы пересмотреть после фоновой проверки; не присваивать по ним статус BGC.
- До адресной загрузки контекстов проверить controls и число hits. Если `≥95.31` дают более 500 уникальных белков, остановиться и откалибровать domain coverage/AsnB-негативы, не ужесточать порог ради удобного ответа.
- Для всех high hits и воспроизводимой выборки grey hits извлекать ДНК/GFF только из соответствующих геномов; проверять A/B/C/(D), край контига и качество MAG. Находки ниже completeness 90% или с contamination ≥5% не исключать из diversity, но не приоритизировать для HGT без дополнительных оснований.
- Один maize-каталог не отвечает на cross-environment routes. После его результата решить, нужен ли один ограниченный внешний поиск, не начинать глобальную выгрузку автоматически.

## Что проверить человеку

1. Согласны ли вы с границей проекта: «профильный C-hit + BGC-контекст», с отдельными canonical/noncanonical ветвями?
2. Считать ли *P. asymbiotica* полезным provisional-кейсом из-за согласованности A–D, несмотря на неустойчивое C-дерево?
3. Принимать ли описанный scope S3 и правило остановки `>500` sequence hits?

## Артефакты

- Нормализация: [`literature_reference_summary.json`](../results/tables/literature_reference_summary.json), [`literature_bgc_manifest.tsv`](../data/processed/literature_bgc_manifest.tsv).
- Профили: [`literature_profile_calibration_summary.json`](../results/tables/literature_profile_calibration_summary.json), [`literature_profile_cv_scores.tsv`](../results/tables/literature_profile_cv_scores.tsv).
- Устойчивость дерева: [`focal_split_sensitivity.tsv`](../results/tables/focal_split_sensitivity.tsv).
- A–D: [`pseudomonadota_core_summary.json`](../results/tables/pseudomonadota_core_summary.json), [`p_asym_core_affinity.tsv`](../results/tables/p_asym_core_affinity.tsv).
- Полные контексты: [`focal_context_manifest.tsv`](../results/tables/focal_context_manifest.tsv), [`focal_context_features.tsv`](../results/tables/focal_context_features.tsv).
- Полная мотивация и история изменений: [`research-log.md`](research-log.md).
