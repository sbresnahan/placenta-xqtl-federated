chr22 phenotyping reference fixture (tests only; rebuilt by tests/fixtures/make_refs_fixture.sh).
annotation/  GENCODE v45 chr22 stand-in for HPLRv2, normalized like the atlas (prep_PANTRY_gtf.R form)
txrevise/    txrevise events from that annotation (PANTRY txrevise scripts, atlas settings)
editing/     REDIportal hg38 sites on chr22 (BED6), cut from the MDACC reference record
music/       mini MuSiC reference: 254 chr22 genes, <= 30 cells per cell type x donor (3,857 cells)
             editing/ and music/ are made by tests/fixtures/extract_atlas_refs_chr22.sh (seadragon)
The genome is not stored here: the tests cut chr22 out of the image's GRCh38 FASTA.
