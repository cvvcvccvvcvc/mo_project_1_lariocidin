# Лариоцидиноподобные BGC: поиск и проверка HGT

Курсовое исследование состоит из двух самостоятельных расчётов. [Diversity](diversity/README.md) ищет и курирует лариоцидиноподобные BGC в фиксированных публичных выборках; [HGT](hgt/README.md) проверяет филогенетические гипотезы на собственном снимке входов. Исследовательский вопрос и границы вывода — в [научной рамке](docs/research.md). Предсказанный BGC не доказывает синтез продукта; текущие данные HGT не устанавливают событие переноса.

## Повторить расчёт

Нужны Linux, micromamba, доступ к публичным источникам и для полного анализа кластер со Slurm. Код и небольшие зафиксированные таблицы находятся в Git; геномы, окружения, журналы и крупные архивы не коммитятся. Запишите commit и не смешивайте его результаты с другой версией выборки.

1. После клонирования перейдите в `diversity/`. Создайте окружение: `micromamba create -y -p .env -f environment/environment.linux-64.lock` (на macOS — `environment/environment.osx-arm64.lock`; общий рецепт — `environment/environment.yml`).
2. Проверочный анализ запускается одной командой: `bash scripts/run_reference.sh`. Затем выполните `.env/bin/python -m unittest discover -s tests -v`. Этот маршрут проверяет установки, входы и контрольный reference; он не заменяет полный поиск.
3. Для полного поиска получите два архива, перечисленных в [реестре внешних входов](diversity/inputs/provenance/off_git_assets.tsv), поместите их в `diversity/inputs/archives/` и из `diversity/` проверьте `sha256sum -c inputs/provenance/off_git_assets.sha256`. Архивы пока не опубликованы: одного Git clone для исторического полного повтора недостаточно.
4. На кластере из `diversity/` последовательно отправьте `bash scripts/submit_baseline.sh`, затем после успешного завершения задания `bash scripts/submit_species.sh`, `bash scripts/submit_additional.sh extension`, `bash scripts/submit_additional.sh controls` и после готовности всех 17 shards `bash scripts/submit_diversity_review.sh`. Эти входные команды задают зависимости Slurm. Подробная последовательность, входы и проверяемые выходы — в [методе diversity](diversity/docs/method.md).
5. Сравните итоговые таблицы с [зафиксированным результатом](diversity/docs/results.md) и оформите новый manifest с commit, входами, версиями и SHA-256. [HGT](hgt/README.md) запускается отдельно на своём снимке; нельзя молча подменять его новым diversity-каталогом.

## Кластер без личных данных в Git

Оператор получает адрес login node, логин, ключ, при необходимости jump host и рабочий каталог у администратора. Их хранят вне репозитория: личный alias `mo-cluster` в `~/.ssh/config`, приватный ключ в `~/.ssh/`, путь к клону в переменной `LARIOCIDIN_WORKDIR`. Типовая сессия после настройки:

```sh
ssh mo-cluster
cd "$LARIOCIDIN_WORKDIR/diversity"
squeue -u "$USER"
bash scripts/submit_baseline.sh
```

На login node готовьте входы и отправляйте задания; тяжёлые расчёты выполняйте через Slurm. Планируйте в пределах 64 CPU, 512 GiB RAM и 250 GB диска суммарно; актуальные partition, QOS и квоту проверяйте на самом кластере. Не добавляйте в Git hostname, имена учётных записей, ключи, личные абсолютные пути или большие данные. Крупные входы сверяются по SHA-256, а не по имени файла.
