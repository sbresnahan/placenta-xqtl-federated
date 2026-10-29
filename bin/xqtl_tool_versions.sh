#!/usr/bin/env bash
# Print the version of every tool and library the pipeline uses, as a TSV
# (tool <TAB> version-or-MISSING). Used by --mode smoke (S02) and recorded in each
# export's provenance. Exit status is always 0; callers inspect the MISSING rows.
set -u
printf "tool\tversion\n"
v() {  # name, command producing the version (first line kept)
    local name=$1; shift
    local out
    if out=$( (eval "$@") 2>&1 | grep -v '^\s*$' | head -1 ) && [ -n "$out" ] && ! echo "$out" | grep -qiE 'command not found|No such file|not installed|there is no package|ModuleNotFoundError|Error'; then
        printf "%s\t%s\n" "$name" "$(echo "$out" | tr '\t' ' ')"
    else
        printf "%s\tMISSING\n" "$name"
    fi
}
v STAR          "STAR --version"
v salmon        "salmon --version"
v samtools      "samtools --version | head -1"
v bcftools      "bcftools --version | head -1"
v htslib_bgzip  "bgzip --version | head -1"
v tabix         "tabix --version | head -1"
v featureCounts "featureCounts -v 2>&1 | grep -i featurecounts"
v gffread       "gffread --version"
v regtools      "regtools --version 2>&1 | grep -i version"
v bedtools      "bedtools --version"
v plink2        "plink2 --version"
v CrossMap      "CrossMap --version 2>&1 | grep -i crossmap"
v IRFinder-S    "IRFinder --version 2>&1 | grep -i 'irfinder\|version'"
v python        "python3 --version"
for d in numpy pandas scipy pyarrow torch tensorqtl pandas-plink pgenlib pysam jsonschema PyYAML; do
    v "py:$d" "python3 -W ignore -c 'import importlib.metadata as m; print(m.version(\"$d\"))'"
done
v R             "R --version | head -1"
for p in edgeR DESeq2 tximport rtracklayer qvalue sva Rhcpp MuSiC GENESIS SNPRelate SeqArray coloc susieR colocboost glmnet data.table optparse; do
    v "R:$p" "Rscript -e 'cat(as.character(packageVersion(\"$p\")))'"
done
exit 0
