# ISO_DAF

Reusable Snakemake workflow for DAF-seq analysis, including QC, alignment,
signal extraction, FiberHMM integration, delivery tables, and offline reports.

The public tree contains reusable rules, environments, tests, and a demo
profile. Real sample sheets, references, customer delivery material, and
result trees are excluded.

~~~bash
make lint
make smoke
make dry-run CONFIG=config/profiles/demo.yaml
~~~

Use a local profile for real data and keep runtime artifacts outside Git.
