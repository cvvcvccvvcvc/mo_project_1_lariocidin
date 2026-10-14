# Скрипты diversity

Запускайте сценарии из [README проекта](../README.md), а не файлы по алфавиту.

| Входная команда | Что делает |
|---|---|
| `submit_baseline.sh` | Отправляет baseline и подготовку diversity через Slurm после проверки архивов. |
| `run_reference.sh` | Небольшой reference-повтор и контроль установки. |
| `fetch_research_inputs.py --diversity` | Восстанавливает зафиксированные metadata snapshots. |
| `submit_species.sh` | Отправляет восемь видовых shards, контексты и сборку с зависимостями. |
| `submit_additional.sh extension|controls` | Отправляет межсредовой этап или контроли со своими выборками. |
| `submit_diversity_review.sh` | После готовых 17 shards ставит повторную курацию и сборку через Slurm. |
| `export_diversity_handoff.py --output <каталог>` | Делает новый версионированный снимок для HGT. |

`prepare_*` и `fetch_*` готовят источники, `screen_*`/`extract_*` выполняют поиск и контексты, `collect_*`/`merge_*`/`compare_*` сводят каталоги, `audit_*`/`validate_*` проверяют основания. `slurm/` — задания массивов; `.sbatch` и shard-скрипты вызываются из конкретного этапа. Скрипты раннего reference/Maize этапа оставлены как зависимости воспроизведения. Устаревшие расчёты S5 HGT удалены: текущая HGT-линия живёт в `../../hgt/`.
