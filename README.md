# HPDB data branch

This orphan branch hosts a compact, browser-friendly build of the Hypothetical
Protein Database (HPDB) so it can be loaded over a CORS-enabled raw URL. It is
kept separate from the main source branch to avoid storing the binary in the
code history.

The file on this branch is a lite build intended for in-browser exploration, not
the full database. To keep it under the GitHub file size limit, the verbose
evidence and rationale columns are omitted, only the strongest evidence rows per
protein and channel are retained, functional terms are capped per protein and
term type, and only the nearest genomic neighbours are kept. The core protein,
family and structure tables are retained in full. The complete database with all
evidence rows and full metadata is published as a GitHub release asset.

Open the lite database in your browser with Datasette Lite:

https://lite.datasette.io/?url=https%3A%2F%2Fraw.githubusercontent.com%2FLLQ95%2FPUFprobe%2Fhpdb-data%2Fhpdb.sqlite

Download the full database from the release page:

https://github.com/LLQ95/PUFprobe/releases
