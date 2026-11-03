# HPDB data branch

This orphan branch hosts the built Hypothetical Protein Database (HPDB) SQLite
file so it can be loaded over a CORS-enabled raw URL. It is intentionally kept
separate from the main source branch to avoid storing the large binary in the
code history.

Open the database in your browser with Datasette Lite:

https://lite.datasette.io/?url=https%3A%2F%2Fraw.githubusercontent.com%2FLLQ95%2FPUFprobe%2Fhpdb-data%2Fhpdb.sqlite
