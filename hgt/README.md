# Отдельная HGT-линия

Для запуска на кластере сначала выполните [общую подготовку и проверку доступа](../README.md#доступ-и-запуск-в-кластере). HGT ведётся в отдельном checkout и рабочей области. Перед `sbatch` экспортируйте `HGT_ROOT` (каталог этой линии), `HGT_BIN` (каталог `bin` зафиксированного окружения), а для H2 также `HGT_GENOME_ROOT` (замороженные H1-геномы). Эти личные пути задаются оператором вне Git; сценарии прекращают работу, если необходимой переменной нет. Нужны исходные файлы из манифестов `config/` с проверенными SHA-256. Большие `inputs/` и промежуточные `results/` в Git не хранятся, поэтому полный HGT-повтор требует отдельной публикации этих входов.

Текущий дизайн и состояние: [план](plan.md). Научное основание выбора методов и выполненные проверки: [обзор 13.10.2026](reviews/2026-10-13-foundation-and-methods.md). Решение пользователя: canonical и noncanonical исследуются как две самостоятельные линии.

Все новые записи принадлежат этой директории. Общие документы, diversity-код, входы и результаты другой линии не изменялись. Снимки исходных архивов и прежнего S5 находятся в `inputs/2026-10-13/`; их происхождение и SHA-256 указаны в `config/input_snapshot_2026-10-13.json`. Снимок общих документов отражает состояние на момент копирования и не заменяет их активные версии.

Воспроизвести выполненную инвентаризацию из корня репозитория:

```sh
projects/lariocidin/.env/bin/python projects/lariocidin/hgt/scripts/audit_published_inputs.py
bash projects/lariocidin/hgt/scripts/run_waafle_fixture.sh
```

Для первого скрипта нужны Python и Biopython; для WAAFLE-пилота дополнительно NumPy. Использовано уже существовавшее окружение, без установки в него новых пакетов. Его версии сохранены в `config/preparation_environment.json`. Скрипту WAAFLE можно передать другой совместимый Python первым аргументом.

WAAFLE-пилот воспроизводит gene calling и scoring авторского демо с готовыми BLAST hits. Координаты генов совпадают с авторским GFF. Scoring **не воспроизводит** сохранённые авторские таблицы; расхождение описано в обзоре, оно не объявлено успешным тестом равенства. Ни поиск по большой базе, ни проверка reads, ни анализ наших BGC здесь не выполнялись.

Завершённые находки diversity-линии позднее принимаются отдельным новым снимком по контракту в плане. Текущие снимки не перезаписываются. H2 выполняется: component trees, три Bacteria74 host trees и отдельная Nature-основа завершены и сопоставлены. [Разбор готовых результатов обеих линий](reviews/2026-10-14-H2-ready-tree-comparison.md) содержит шесть вручную выбранных отношений с независимой gCF и ограничениями. Основное Actino138 host tree продолжает1000UFBoot+BNNI на64CPU; полное сравнение host-панелей ожидается. DTL, AU и каталог подтверждённых HGT-событий ещё не выполнены.

H1-аудит опубликованных локусов: [результаты и ограничения](reviews/2026-10-13-H1-carriers-and-relatives.md), [каталог носителей](results/tables/H1_carrier_catalogue.tsv), [компоненты-кандидаты](results/tables/H1_component_candidates.tsv), [сборки](results/tables/H1_genome_quality.tsv), [реализованная панель родственников](results/tables/H1_relative_panel.tsv). Проверены 320 геномных входов; 251/342 экспорта связаны с носителем. В исходном H1 четыре Nature-системы сохранены как неопределённые; H2 подтвердил их геномную noncanonical принадлежность по precursor+C, сохранив отдельную диагностику и неизвестную биохимическую функцию. Первичная связь с последовательностью не отменяет последующую ортологическую курацию и marker QC.

Подготовка из фиксированных архивов и метаданных:

```sh
projects/lariocidin/.env/bin/python projects/lariocidin/hgt/scripts/prepare_h1.py
bash projects/lariocidin/hgt/scripts/run_h1_homology.sh
projects/lariocidin/.env/bin/python projects/lariocidin/hgt/scripts/fetch_h1.py
projects/lariocidin/.env/bin/python projects/lariocidin/hgt/scripts/prepare_h1_supplement.py --download
```

`fetch_h1.py` использует сохранённые Eutils-ответы и загружает только отсутствующие файлы; изменения URL/содержимого уже зарегистрированного снимка приводят к ошибке. NPDC-кандидаты теперь получаются тем же подготовительным скриптом. Первичная таблица подбора родственников воспроизводится из неизменного `inputs/H1_2026-10-13/relative_selection_basis/`, а не из текущего расширенного каталога: это исключает обратную зависимость выборки от результатов. Первоначальные metadata/source snapshots нужны для полного повторения; их SHA-256 и происхождение сохранены в manifests. Выбранные полные геномы имеют версионированные NCBI URL.

Основной аудит следует запускать через Slurm в отдельной HGT-директории, по `scripts/H1_correspondence.sbatch` и `scripts/H1_supplement.sbatch`; B1/B2-поиск — `scripts/H1_B_homology.sbatch`, восстановление BGSC 4BL1 — `scripts/H1_4BL1.sbatch`. Последний использует отдельные зафиксированные input table/metadata из NCBI поиска штамма. Сценарии не запускают расчёты на login node и не устанавливают пакеты в общее окружение. Перед повторной отправкой переносить только собственные файлы согласно `config/H1_*_input_manifest.json` и проверять соответствующие `.sha256`. Paths/runtime указаны для разрешённого ITMO host; при другом месте расчёта заменить их явно.

Выполненные задания: `937682`, `937712`, `937720`, `937740`; accounting и manifests находятся в `config/`. Все завершились `0:0`. Первый запуск использовал сохранённую в `inputs/H1_2026-10-13/run_937682/` версию аудитора; текущая версия добавляет только выбор отдельной input table/metadata для дополнений. Исходные версии изменившегося NCBI manifest сохранены в `run_937682/` и `run_937712/`. Данные и результаты старого S5 остаются отдельно.

Сборка итогового каталога после возврата завершённых shards:

```sh
projects/lariocidin/.env/bin/python projects/lariocidin/hgt/scripts/collect_h1.py \
  --shard-root projects/lariocidin/hgt/results/H1_shards \
  --shard-root projects/lariocidin/hgt/results/H1_supplement_shards \
  --shard-root projects/lariocidin/hgt/results/H1_4BL1_shards
```

Нужны Biopython и BLAST+; основные версии Python/Biopython/NumPy совпадают с `config/preparation_environment.json`, BLASTP на кластере — 2.17.0+. Генные назначения и исключения читаются из `config/H1_lineage_review.tsv` и `H1_named_anchor_review.tsv`; первичные последовательности при курации не меняются. Контроли пока не считаются BGC-negative, сборки не прошли проверку reads/контаминации, а HGT ещё не оценён.

H2: [методы и ортология](reviews/2026-10-13-H2-orthology-and-phylogeny.md), [готовые результаты](reviews/2026-10-14-H2-ready-tree-comparison.md). Основной Actino138 host fit ещё не завершён. Выполнены 74-marker QC и дополнительная авторская Actinobacteria138 панель, проверка C и контекстных компонентов, отдельная Nature-курация. [Выбор моделей](config/H2_host_model_decision.json) утверждён пользователем; [соответствие финальным выравниваниям](config/H2_final_alignment_model_correspondence.json) проверено через Slurm. Незавершённый MFP+MERGE сохранён с checkpoints/logs как ресурсный пилот, без утверждения равнозначности схем. Все вычисления используют HGT snapshots и собственные результаты.

`scripts/submit_H2_collect.sh <marker-search-job-id>` ставит сборку H2 только после успешного завершения переданного Slurm job; ID берите из своего текущего запуска, исторические номера заданий не используются.

Для повторения H2 переносить только HGT files согласно frozen input manifests, задавать `HGT_ROOT` на отдельный checkout и `HGT_GENOME_ROOT` на фиксированные H1 genomes. Read-only runtime задаётся `HGT_BIN`; Prodigal/trimAl — собственным `HGT_SMALL_BIN` по `config/H2_tools_explicit.txt`. Основные этапы через Slurm: `prepare_h2.py` → `H2_native.sbatch` → `H2_marker_collect.sbatch`/`H2_components.sbatch` → отдельные marker jobs → `h2_prepare_trees.py --hosts --fixed-host-models` → `H2_host_tree.sbatch` с `HGT_FIXED_HOST_MODELS=1`. Для дополнительной панели используются `H2_actino_search.sbatch` и `--actino`/`HGT_ACTINO=1`; отдельные host/component arrays индексируются зафиксированными `results/H2/tree_jobs.tsv` и `results/H2/actino/tree_jobs.tsv`. Не переписывать порядок этих таблиц во время active array runs. Проверка финальных моделей — `H2_verify_final_models.sbatch`; итоговые gCF и сопоставления/рисунки — `H2_concordance.sbatch`, `H2_compare_plot.sbatch` только после успешного завершения основных trees. Individual marker/domain/extra-copy diagnostics с `-fast` не являются полными основными ML trees.

Для разбора уже завершённых наборов `HGT_READY_ONLY=1` в `H2_compare_plot.sbatch` пишет отдельный `comparisons_partial` и пропускает fits без итогового task summary. Он проверяет leaf IDs и1000bootstrap trees для completed full ML, создаёт описательные сравнения, извлекает свидетельства по вручную зафиксированным `H2_review_relation_queries.json` и рисует trees. `H2_concordance.sbatch` использует topology-only unrooted marker copies; исходные lengths/supports не меняет. Дополнительный небольшой `h2_extra_C_context.py` сравнивает четыре предварительно выбранных C±5kb участка BLASTN, проверяет нативные CDS и точные опубликованные precursor alleles; не классифицирует функцию и не меняет основные входы. Его fixed-input outputs находятся в `results/H2/orthology_diagnostic/context_check/`.
