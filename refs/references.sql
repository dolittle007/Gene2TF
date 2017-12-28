-- NCBI refseq
CREATE TABLE refseq (ID INTEGER PRIMARY KEY AUTOINCREMENT,chr VARCHAR(10) NOT NULL,start INT NOT NULL, end INT NOT NULL,transcript VARCHAR(100) NOT NULL, symbol VARCHAR(100) NOT NULL, strand CHAR(1), geneid VARCHAR(100) NOT NULL, biotype VARCHAR(100) NOT NULL);
CREATE TABLE tmp (chr VARCHAR(10) NOT NULL,start INT NOT NULL,end INT NOT NULL,transcript VARCHAR(100) NOT NULL, symbol VARCHAR(100) NOT NULL,strand CHAR(1),geneid VARCHAR(100) NOT NULL,biotype VARCHAR(100) NOT NULL);
.separator "\\t"
.mode tabs
.import Homo_sapiens.GRCh37.refseq.bed tmp
INSERT INTO refseq (chr, start, end, transcript, symbol, strand, geneid, biotype) select chr, start, end, transcript, symbol, strand, geneid, biotype from tmp;
DROP TABLE tmp;
.header on
CREATE INDEX transcript_index_refseq ON refseq (transcript);
CREATE INDEX symbol_index_refseq ON refseq (symbol);
CREATE INDEX geneid_index_refseq ON refseq (geneid);
CREATE INDEX biotype_index_refseq ON refseq (biotype);
-- EBI ensembl
CREATE TABLE ensembl (ID INTEGER PRIMARY KEY AUTOINCREMENT,chr VARCHAR(10) NOT NULL,start INT NOT NULL, end INT NOT NULL,transcript VARCHAR(100) NOT NULL, symbol VARCHAR(100) NOT NULL, strand CHAR(1), geneid VARCHAR(100) NOT NULL, biotype VARCHAR(100) NOT NULL);
CREATE TABLE tmp (chr VARCHAR(10) NOT NULL,start INT NOT NULL,end INT NOT NULL,transcript VARCHAR(100) NOT NULL, symbol VARCHAR(100) NOT NULL,strand CHAR(1),geneid VARCHAR(100) NOT NULL,biotype VARCHAR(100) NOT NULL);
.separator "\\t"
.mode tabs
.import Homo_sapiens.GRCh37.ensembl.bed tmp
INSERT INTO ensembl (chr, start, end, transcript, symbol, strand, geneid, biotype) select chr, start, end, transcript, symbol, strand, geneid, biotype from tmp;
DROP TABLE tmp;
.header on
CREATE INDEX transcript_index_ensembl ON ensembl (transcript);
CREATE INDEX symbol_index_ensembl ON ensembl (symbol);
CREATE INDEX geneid_index_ensembl ON ensembl (geneid);
CREATE INDEX biotype_index_ensembl ON ensembl (biotype);
