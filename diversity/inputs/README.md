# Входы и provenance

`selections/` содержит общие списки источников и решения по выборке; `selections/shards/` — 17 зафиксированных частей видового, межсредового и контрольного поиска. `provenance/` содержит версии, URL, run manifests и SHA-256; `provenance/shards/` — 34 журнала загрузок по тем же частям. Номер shard относится к выборке, а не к порядку запуска скриптов.

`archives/` хранит исходные снимки: крупные файлы исключены из Git и перечислены в `provenance/off_git_assets.tsv`; перед полным повтором проверьте `provenance/off_git_assets.sha256`. `provenance/historical/` сохраняет исходные manifests прошлых запусков с путями и хешами той версии кода; для нового запуска создавайте новый manifest.

Префиксы `diversity_species_*`, `diversity_extension_*` и `diversity_controls_*` отделяют три этапа. Старые `maize_*`, `s3_*`, `s5_*` и `reference_*` относятся к ограниченному baseline и не заменяют текущие выборки.
